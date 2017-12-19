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
FILTER_SHAPE = [3, 3, 1, 2]
FILTER = np.zeros(shape=(3, 3, 1, 2), dtype=np.float32) # 2 features
FILTER[:,:,0,0] = 1/6
FILTER[:,:,0,1] = 1
FILTER[0,0,0,:] = 0.
FILTER[-1,-1,0,:] = 0.

CONV_STRIDES = [1,1,1,1]
POOL_KSIZE = [1,2,2,1]
POOL_STRIDES = [1,2,2,1]
FULLY_CONNECTED_SIZE = []

BATCH_SIZE = 50
IM_WIDTH = 40
IM_HEIGHT = 39
LEARNING_RATE = 0.005
EVAL_SIZE = 200
GENERATIONS = 500
CONV_FEATURES = 25


tf.reset_default_graph()
sess = tf.Session()

##########################
###### READ THE FILE #####
##########################
def PRRevaluation(**kwargs):
    
    f_images = kwargs['images']
    f_labels = kwargs['labels']
    
    #############################################
    ####### EVALUATION and RETURN TENSORS #######
    #############################################
    with open(f_images, 'rb') as f:
        images = pic.load(f)
    with open(f_labels, 'rb') as ff:
        labels = pic.load(ff)

    # devide in to train and test
    images_train = images[:int(len(images)*0.6)]
    labels_train = labes[:int(len(images)*0.6)]
    images_test = images[int(len(images)*0.6):]
    labels_test = labels[int(len(images)*0.6):]

    # define placeholdes fot train and test
    train_im_input_shape = (BATCH_SIZE, IM_WIDTH, IM_HEIGHT, 1)
    train_im_input = tf.placeholder(tf.float32, shape=train_im_input_shape)
    train_target = tf.placeholder(tf.float32, shape=(BATCH_SIZE))

    test_im_input_shape = (EVAL_SIZE, IM_WIDTH, IM_HEIGHT, 1 )
    test_im_input = tf.placeholder(tf.float32, shape=test_im_input_shape)
    test_target = tf.placeholder(tf.float32, shape=(EVAL_SIZE))

    # Convolutional layer variables
    conv_weight = tf.Variable(tf.truncated_normal(FILTER_SHAPE, stddev=0.1,
                                                   dtype=tf.float32))
    conv_bias = tf.Variable(tf.zeros([FILTER_SHAPE[-1]], dtype=tf.float32))


    #### test
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
