#***********************************************************************
# Copyright (C) 2017 the Imaging X-ray Polarimetry Explorer (IXPE) team.
#
# For the license terms see the file LICENSE, distributed along with this
# software.
# 
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#***********************************************************************


import os
import imp
import argparse
import numpy as np
import pickle as pic
import matplotlib.pyplot as plt
import tensorflow as tf

from matplotlib import cm, colors
cmap = cm.get_cmap('viridis')


__description__ = 'Data Transformation: from hexagonal to squared pixels to tensors'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-im', '--images', type=str, required=True,
                    help='images file')
PARSER.add_argument('-lb', '--labels', type=str, required=True,
                    default=None, help='labels file')


###################################
####### PATTERN-REC PARAM #########
###################################
FILTER_W = 3
FILTER_H = 3
FEATURES = 1
FILTER_SHAPE = [FILTER_W, FILTER_H, 1, FEATURES]


CONV_STRIDES = [1,1,1,1]
POOL_KSIZE = [1,2,2,1]
POOL_STRIDES = [1,2,2,1]
FULLY_CONNECTED_SIZE = []

BATCH_SIZE = 50
IM_WIDTH = 40
IM_HEIGHT = 39
LEARNING_RATE = 1e-4
EVAL_SIZE = 200
EVAL_EVERY = 50
GENERATIONS = 1000
###################################

def conv2d(x, W, padding):
    return tf.nn.conv2d(x, W, strides=CONV_STRIDES, padding=padding)

def max_pool_2x2(x):
    return tf.nn.max_pool(x, ksize=POOL_KSIZE, strides=POOL_STRIDES, padding='SAME')

def weight_variable(shape):
    #initial = tf.truncated_normal(shape, mean=0.0, stddev=1.)
    #initial = tf.random_uniform(shape, minval=-1, maxval=1)
    initial = tf.random_normal(shape)
    return tf.Variable(initial)

def bias_variable(shape):
    initial = tf.constant(0.01, shape=shape)
    return tf.Variable(initial)

def next_data_batch(a, b, index, batchsize):
    return np.array([a[index-batchsize:index]]), np.array(b[index-batchsize:index])

def simultaneous_shuffle(a, b):
    shuffled_a = np.empty(a.shape, dtype=a.dtype) # create an empty array
    shuffled_b = np.empty(b.shape, dtype=b.dtype) # create an empty array
    permutation = np.random.permutation(len(a)) # create a permutation
    for old_index, new_index in enumerate(permutation): # loop over the array
        shuffled_a[new_index] = a[old_index] # create the permutation
        shuffled_b[new_index] = b[old_index] # create the permutation
    return shuffled_a, shuffled_b # return

def nearest_neighbor_kernel(shape):
    a_list = []
    b_list = []
    weights = tf.zeros(shape)
    b_list.append(weights)
    weights = weight_variable(shape)
    b_list.append(weights)
    weights = weight_variable(shape)
    b_list.append(weights)
    a_list.append(tf.stack(b_list)) # pack the first row together
    b_list = []
    weights = weight_variable(shape)
    b_list.append(weights)
    weights = weight_variable(shape)
    b_list.append(weights)
    weights = weight_variable(shape)
    b_list.append(weights)
    a_list.append(tf.stack(b_list)) # pack the first and second row together
    b_list = []
    weights = weight_variable(shape) # initialize weights with shape (input_dimension, output_dimension)
    b_list.append(weights)
    weights = weight_variable(shape)
    b_list.append(weights)
    weights = tf.zeros(shape) # does not belong to hexagon; generate zeros
    b_list.append(weights)
    a_list.append(tf.stack(b_list)) # pack the first, second and third row together
    hexKernel = tf.stack(a_list)
    return hexKernel

def my_conv_net(input_data):
    W_conv1 = nearest_neighbor_kernel([1,FEATURES])
    b_conv1 = bias_variable([FEATURES])
    h_conv1 = tf.nn.softsign(tf.add(conv2d(input_data, W_conv1, "SAME"), b_conv1))
    h_pool1 = max_pool_2x2(h_conv1)
    # print(h_pool1.shape) # [?, 20, 20, FEATURES]
    W_conv2 = nearest_neighbor_kernel([FEATURES,FEATURES])
    b_conv2 = bias_variable([FEATURES])
    h_conv2 = tf.nn.softsign(tf.add(conv2d(h_pool1, W_conv2, "SAME"),b_conv2))
    h_pool2 = max_pool_2x2(h_conv2)
    # print(h_pool2.shape) # [?, 10, 10, FEATURES]
    W_conv3 = nearest_neighbor_kernel([FEATURES,FEATURES])
    b_conv3 = bias_variable([FEATURES])
    h_conv3 = tf.nn.softsign(tf.add(conv2d(h_pool2, W_conv3, "SAME"), b_conv3))
    h_pool3 = max_pool_2x2(h_conv3)
    # print(h_pool3.shape) # [?, 5, 5, FEATURES]
    W_fc1 = weight_variable([5 * 5 * FEATURES, 1])
    b_fc1 = bias_variable([1])
    h_pool3_flat = tf.reshape(h_pool3, [-1, 5 * 5 * FEATURES])
    h_fc1 = tf.nn.softsign(tf.add(tf.matmul(h_pool3_flat, W_fc1), b_fc1))
    W_fc2 = weight_variable([1, 1])
    b_fc2 = bias_variable([1])
    model = tf.matmul(h_fc1, W_fc2) + b_fc2
    return model

def get_label_normalized(labels_arr):
    N = len(labels_arr)
    mu = np.mean(labels_arr)
    sigma = np.sqrt(1 / (N - 1) * np.sum(((labels_arr - mu)**2)))
    labels_norm = (labels_arr - mu) / sigma
    return labels_norm

def get_image_normalized(images_arr):
    N = len(images_arr)
    mu_ = np.mean(images_arr, axis=0)
    sigma_ = np.sqrt(1 / (N - 1) * np.sum(np.subtract(images_arr,mu_)**2, axis=0))
    images_norm = np.divide(np.subtract(images_arr, mu_)[0], sigma_[0])
    return images_norm

##########################
###### READ THE FILE #####
##########################
def PRRevaluation(**kwargs):
    
    tf.reset_default_graph()
    sess = tf.Session()
    
    f_images = kwargs['images']
    f_labels = kwargs['labels']
    
    #############################################
    ####### EVALUATION and RETURN TENSORS #######
    #############################################
    with open(f_images, 'rb') as f:
        images = pic.load(f)
    with open(f_labels, 'rb') as ff:
        labels = pic.load(ff)

    index_clean = np.where(labels[:-2,0] != -1000.)[0]
    images = images[index_clean]
    labels = labels[index_clean]

    height, width = images[0].shape

    # devide in to train and test
    images_train = images[:int(len(images)*0.7)]
    labels_train = labels[:int(len(images)*0.7)][:,1]
    images_test = images[int(len(images)*0.7):]
    labels_test = labels[int(len(images)*0.7):][:,1]

    labels_train_norm = get_label_normalized(labels_train)

    x_image = tf.placeholder(tf.float32, shape=[None, height, width, 1])
    y_ = tf.placeholder(tf.float32, shape=[None, 1])

    eval_x_image = tf.placeholder(tf.float32, shape=[None, height, width, 1])
    eval_y_ = tf.placeholder(tf.float32, shape=[None, 1])
    # print(x_image.shape) # [?, 40, 39, 1]

    model_output = my_conv_net(x_image)
    test_model_output = my_conv_net(eval_x_image)

    loss = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(logits=model_output, labels=y_))


    train_step = tf.train.AdamOptimizer(LEARNING_RATE).minimize(loss)

    init = tf.global_variables_initializer()
    sess.run(init)

    train_loss = []

    for i in range(GENERATIONS):
        rand_index = np.random.choice(len(images_train), size=BATCH_SIZE)
        train_phis, train_pixels = np.transpose([labels_train_norm[rand_index]]), images_train[rand_index]
        train_pixels_expanded = tf.expand_dims(train_pixels, -1).eval(session=sess)
        train_dict = {x_image: train_pixels_expanded, y_: train_phis}

        sess.run(train_step, feed_dict=train_dict)
        temp_train_loss = sess.run(loss, feed_dict=train_dict)
    
        if (i+1) % EVAL_EVERY == 0:
            eval_index = np.random.choice(len(images_test), size=EVAL_SIZE)
            eval_x = images_test[eval_index]
            eval_x = np.expand_dims(eval_x, 3)
            eval_y = np.transpose([labels_test[eval_index]])
            test_dict = {eval_x_image: eval_x, eval_y_: eval_y}
            test_preds = sess.run(test_model_output, feed_dict=test_dict)
            
            # Record and print results
            train_loss.append(temp_train_loss)
            acc_and_loss = [(i+1), temp_train_loss]
            print('Generation # {}. Train Loss: {:.2f}'.format(*acc_and_loss))

    test_pixels_expanded = tf.expand_dims(images_test, -1).eval(session=sess)
    predicted_phis = sess.run(test_model_output, feed_dict={eval_x_image: test_pixels_expanded, eval_y_:np.transpose([labels_test])})
    #for i in range(0,100):
    #print(labels_test[i], predicted_phis[i])

    plt.hist(labels_test, bins=20)
    plt.hist(predicted_phis, bins=20)
    plt.show



    #### some plotting
    
    FILTER = np.zeros(shape=(3, 3, 1, 2), dtype=np.float32) # 2 features
    FILTER[:,:,0,0] = 1/6
    FILTER[:,:,0,1] = 1
    FILTER[0,0,0,:] = 0.
    FILTER[-1,-1,0,:] = 0.

    events = len(images)
    height, width = images[0].shape
    images = images.reshape(events, height, width, 1)
    print('input shape:', images.shape)

    Xc = tf.placeholder(tf.float32,
                       shape=(None, height, width, 1))
    feature_maps = tf.constant(FILTER)
    convolution = tf.nn.conv2d(
                           Xc,
                           feature_maps,
                           strides=CONV_STRIDES,
                           padding="SAME", 
                           use_cudnn_on_gpu=False)

    output_conv = convolution.eval(session=sess, feed_dict={Xc: images})
    print('convolution output shape:', output_conv.shape)

    Xp = tf.placeholder(tf.float32,
                       shape=output_conv.shape)
    max_pool = tf.nn.max_pool(
                          Xp,
                          ksize=POOL_KSIZE,
                          strides=POOL_STRIDES,
                          padding="VALID")
    output_maxpool = sess.run(max_pool, feed_dict={Xp: output_conv})
    print('maxpooling output shape:', output_maxpool.shape)


    #########################
    ###### SOME DRAWING #####
    #########################
    f, (ax1, ax2) = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 3))
    plt.title('Filter features')
    ax1.imshow(FILTER[:, :, 0, 0], vmin=0, vmax=1)
    ax2.imshow(FILTER[:, :, 0, 1], vmin=0, vmax=1)
    ax1.axis('off')
    ax2.axis('off')

    fig = plt.figure(figsize=(10, 3))
    grid = plt.GridSpec(1, 3) # hspace=0.8, wspace=1.2
    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1])
    ax3 = fig.add_subplot(grid[2])

    ax1.imshow(images[1,:,:,0], cmap=cmap)
    ax1.set_title('MC energy = %.2f KeV, MC phi = %.2f rad'%(labels[1][0], labels[1][1]), size=8)
    ax1.set_aspect('equal')
    ax1.set_anchor('C')
    ax1.axis('off')

    ax2.imshow(output_conv[1,:,:,0], cmap=cmap)
    ax2.set_title('After Convolution layer 1, filter 1', size=8)
    ax2.set_aspect('equal')
    ax2.set_anchor('C')
    ax2.axis('off')

    ax3.imshow(output_maxpool[1,:,:,0], cmap=cmap)
    ax3.set_title('After Pooling layer 1', size=8)
    ax3.set_aspect('equal')
    ax3.set_anchor('C')
    ax3.axis('off')


    plt.show()


if __name__ == '__main__':
    args = PARSER.parse_args()
    PRRevaluation(**args.__dict__)
