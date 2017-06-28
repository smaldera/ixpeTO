import tensorflow as tf
import numpy as np
import h5py as h5

def simultaneous_shuffle(a, b):
    shuffled_a = np.empty(a.shape, dtype=a.dtype) # create an empty array
    shuffled_b = np.empty(b.shape, dtype=b.dtype) # create an empty array
    permutation = np.random.permutation(len(a)) # create a permutation
    for old_index, new_index in enumerate(permutation): # loop over the array
        shuffled_a[new_index] = a[old_index] # create the permutation
        shuffled_b[new_index] = b[old_index] # create the permutation
    return shuffled_a, shuffled_b # return

def read_hdf5_dataset(filename):
    h = h5.File(filename, 'r')
    data = h['data']
    return data    

def next_data_batch(a, b, index, batchsize):
    return np.array([a[index-batchsize:index]]), np.array(b[index-batchsize:index])

def nearest_neighbor_kernel(shape):
        
        a_list = []
        b_list = []
    
        weights = weight_variable(shape) # initialize weights with shape (input_dimension, output_dimension)
        b_list.append(weights)
        weights = weight_variable(shape)
        b_list.append(weights)
        weights = tf.zeros(shape) # does not belong to hexagon; generate zeros
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
    
        weights = tf.zeros(shape)    
        b_list.append(weights)
        weights = weight_variable(shape)
        b_list.append(weights)
        weights = weight_variable(shape)
        b_list.append(weights)
    
        a_list.append(tf.stack(b_list)) # pack the first, second and third row together
    
        hexKernel = tf.stack(a_list) 
        
        return hexKernel # return the kernel

labels = read_hdf5_dataset("labels.hdf5")
pixels = read_hdf5_dataset("pixels.hdf5")
phis = labels[:,5]

train_pixels = pixels[0:10000]
test_pixels = pixels[10000:]

train_phis = phis[0:10000]
test_phis = phis[10000:]

(train_phis.shape[0])


index = 0
batchsize = 64

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.001)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.001, shape=shape)
  return tf.Variable(initial)

def conv2d(x, W, padding):
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding=padding)

def max_pool_2x2(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME')

x_image = tf.placeholder(tf.float32, shape=[None, 40, 40, 1])
y_ = tf.placeholder(tf.float32, shape=[None, 1])

W_conv1 = nearest_neighbor_kernel([1,32])
b_conv1 = bias_variable([32])

h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1, "SAME") + b_conv1)
h_pool1 = max_pool_2x2(h_conv1)

W_conv2 = nearest_neighbor_kernel([32,32])
b_conv2 = bias_variable([32])

h_conv2 = tf.nn.relu(conv2d(h_pool1, W_conv2, "SAME") + b_conv2)
h_pool2 = max_pool_2x2(h_conv2)

W_conv3 = nearest_neighbor_kernel([32,32])
b_conv3 = bias_variable([32])

h_conv3 = tf.nn.relu(conv2d(h_pool2, W_conv3, "SAME") + b_conv3)
h_pool3 = max_pool_2x2(h_conv3)

W_fc1 = weight_variable([5 * 5 * 32, 256])
b_fc1 = bias_variable([256])

h_pool3_flat = tf.reshape(h_pool3, [-1, 5 * 5 * 32])
h_fc1 = tf.nn.relu(tf.matmul(h_pool3_flat, W_fc1) + b_fc1)

W_fc2 = weight_variable([256, 1])
b_fc2 = bias_variable([1])

y_conv = tf.matmul(h_fc1, W_fc2) + b_fc2

loss = tf.nn.l2_loss(y_conv-y_) # define the loss function

train_step = tf.train.AdamOptimizer(1e-4).minimize(loss)

sess = tf.InteractiveSession()

sess.run(tf.global_variables_initializer())

for i in range(10000):
    if (index+batchsize > train_phis.shape[0]):
        index = batchsize
        train_phis, train_pixels = simultaneous_shuffle(train_phis, train_pixels)
    else:
        index += batchsize
    phi_batch, pixel_batch = next_data_batch(train_phis, train_pixels, index, batchsize)
    phi_batch = phi_batch.T
    train_step.run(feed_dict={x_image: pixel_batch, y_: phi_batch})
    print(i, sess.run(loss, feed_dict={x_image: pixel_batch, y_: phi_batch}))

phi_test_batch, pixel_test_batch = next_data_batch(test_phis, test_pixels, 100, 100)
phi_test_batch = phi_test_batch.T
predicted_phis = sess.run(y_conv,feed_dict={x_image: pixel_test_batch, y_: phi_test_batch})
for i in range(0,100):
    print(phi_test_batch[i], predicted_phis[i])

#correct_prediction = tf.equal(tf.argmax(y_conv, 1), tf.argmax(y_, 1))
#accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

#print(accuracy.eval(feed_dict={x: mnist.test.images, y_: mnist.test.labels}))