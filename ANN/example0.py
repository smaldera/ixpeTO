#***********************************************************************
#
# SOFTMAX Model implementation
# (https://www.tensorflow.org/versions/r0.12/tutorials/mnist/beginners/)
#
#************************************************************************


import tensorflow as tf

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_data', one_hot=True)

x = tf.placeholder(tf.float32, shape=[None, 784])
# x isn't a specific value.
# It's a placeholder, a value that we'll input when we ask TensorFlow
# to run a computation. 2-D tensor of floating-point numbers, with a
# shape [None, 784]
# 784 is the number of pixels of the image (flattened)
# None means

W1 = tf.Variable(tf.truncated_normal([784,10], stddev=0.1))
b1 = tf.Variable(tf.zeros([10]))
# A Variable is a modifiable tensor. For machine learning applications,
# one generally has the model parameters be Variables
# 10 is the number of classes (10 numbers)

y = tf.nn.softmax(tf.matmul(x,W1) + b1)
# This implements the model SOFTMAX (y = W * x + b)
# ATT: W and x are flipped in matmul!!!!

y_ = tf.placeholder(tf.float32, [None, 10])
# Typically define what it means for a model to be bad. We call this the
# COST, or the LOSS, it represents how far off our model is from outcome.
# Very common function to determine the loss is the CROSS-ENTROPY (measuring
# how inefficient our predictions are for describing the truth).

cross_entropy = tf.reduce_mean(-tf.reduce_sum(y_ * tf.log(y), reduction_indices=[1]))
# implement the cross-entropy function. tf.reduce_sum adds the elements in
# the second dimension of y, due to the reduction_indices=[1] parameter.
# tf.reduce_mean computes the mean over all the examples in the batch.
"""
    (Note that in the source code, we don't use this formulation, because 
    it is numerically unstable. Instead we apply tf.nn.softmax_cross_entropy_with_logits 
    on the unnormalized logits (e.g., we call softmax_cross_entropy_with_logits 
    on tf.matmul(x, W) + b), because this more numerically stable function
    internally computes the softmax activation. In your code, consider using 
    tf.nn.(sparse_softmax_cross_entropy_with_logits instead)
"""

train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
# TensorFlow knows the entire graph of your computations, it can automatically
# use the backpropagation algorithm to efficiently determine how your variables
# affect the loss you ask it to minimize. Then it can apply your choice of
# optimization algorithm.  The gradient descent algorithm: TensorFlow simply shifts
# each variable a little bit in the direction that reduces the cost.

init = tf.global_variables_initializer()
# operation to initialize the variables we created.

sess = tf.Session()
sess.run(init)
# launch the model in a Session, and now we run the operation that initializes
# the variables

for i in range(1000):
    batch_xs, batch_ys = mnist.train.next_batch(100)
    sess.run(train_step, feed_dict={x: batch_xs, y_: batch_ys})
# Let's train -- we'll run the training step 1000 times!
# Each step of the loop, we get a "batch" of one hundred random data points from
# our training set. We run train_step feeding in the batches data to replace the
# placeholders --> stochastic training

correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
# tf.argmax(y,1) is the label our model thinks is most likely for each input,
# while tf.argmax(y_,1) is the correct label. We can use tf.equal to check if
# our prediction matches the truth. That gives us a list of booleans.

accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
# To determine what fraction are correct, we cast to floating point numbers
# and then take the mean. E.g. [True, False, True, True] would become [1,0,1,1]
# which would become 0.75.

print(sess.run(accuracy, feed_dict={x: mnist.test.images, y_: mnist.test.labels}))
