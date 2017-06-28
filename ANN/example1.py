import tensorflow as tf

#from tqdm import tqdm

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_data', one_hot=True)

x = tf.placeholder(tf.float32, shape=[None, 784])
y_ = tf.placeholder(tf.float32, shape=[None, 10])

W1 = tf.Variable(tf.truncated_normal([784,256], stddev=0.1))
b1 = tf.Variable(tf.zeros([256]))

h1 = tf.nn.relu(tf.matmul(x,W1) + b1)

W2 = tf.Variable(tf.truncated_normal([256,10], stddev=0.1))
b2 = tf.Variable(tf.zeros([10]))

y = tf.matmul(h1,W2) + b2

cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y_, logits=y))

train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)

sess = tf.InteractiveSession()

sess.run(tf.global_variables_initializer())

for _ in range(20000):
  batch = mnist.train.next_batch(100)
  train_step.run(feed_dict={x: batch[0], y_: batch[1]})

correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))

accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

print(accuracy.eval(feed_dict={x: mnist.test.images, y_: mnist.test.labels}))
