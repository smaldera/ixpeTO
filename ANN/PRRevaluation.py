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
import numpy as np
import pickle as pic
import matplotlib.pyplot as plt

from keras.models import Sequential, model_from_json
from keras.layers import Reshape, Dense, Flatten, Conv2D, MaxPooling2D, Dropout, Activation
#from keras.optimizers import RMSprop
#from keras.datasets import mnist
#from keras.utils import np_utils
from keras import initializers, callbacks, regularizers, optimizers, backend

"""
import imp
import argparse

__description__ = 'Data Transformation: from hexagonal to squared pixels to tensors'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-im', '--images', type=str, required=True,
                    help='images file')
PARSER.add_argument('-lb', '--labels', type=str, required=True,
                    default=None, help='labels file')
"""

###################################
####### PATTERN-REC PARAM #########
###################################
BATCH_SIZE = 120
LR = 1e-4
EPOCHS = 200
FRAC = 0.9
SAVE = True
OUTNAME = 'PHI_powerlaw_nopol'
###################################

f_images = '../sim_line59_nopol_images2.pkl'
f_labels = '../sim_line59_nopol_labels2.pkl'

#f_images = '../sim500000_powerlaw_nopol_images2.pkl'
#f_labels = '../sim500000_powerlaw_nopol_labels2.pkl'

#f_images = '../sim_images2.pkl'
#f_labels = '../sim_labels2.pkl'

##################################
class LossHistory(callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
    
    def on_batch_end(self, batch, logs={}):
        self.losses.append(logs.get('loss'))

##########################
###### READ THE FILE #####
##########################
def PRRevaluation():

    with open(f_images, 'rb') as f:
        images_all = pic.load(f)
    with open(f_labels, 'rb') as ff:
        labels_all = pic.load(ff)
    print('Cleaning...')
    index_clean = np.where(labels_all[:,0] != -1000.)[0]
    images = images_all[index_clean]
    labels = labels_all[index_clean][:,-1]
    print('Found %i bad events' %(len(images_all)-len(images)))

    (img_rows, img_cols) = images[0].shape
    input_shape = (1, img_rows, img_cols)
    print('Images shape: %s'%(str(input_shape)))

    perm = np.random.permutation(len(images))
    images = images[perm]
    labels = labels[perm]

    images_train = images[:int(len(images)*FRAC)]
    images_train = images_train.reshape(len(images_train), 1, img_rows, img_cols)
    labels_train = labels[:int(len(images)*FRAC)]

    images_test = images[int(len(images)*FRAC):]
    images_test = images_test.reshape(len(images_test), 1, img_rows, img_cols)
    labels_test = labels[int(len(images)*FRAC):]

    print ('%i events for training'%len(images_train))
    print ('%i events for testing'%len(images_test))

    #############################################
    ####### EVALUATION and RETURN TENSORS #######
    #############################################
    # Convolutional model
    backend.clear_session()
    model = Sequential()
    model.add(Conv2D(64, (3, 3), padding='same', activation='relu',
                 input_shape=input_shape, kernel_regularizer=regularizers.l2(0.01),
                 bias_initializer='RandomUniform', kernel_initializer='RandomUniform'))
    model.add(MaxPooling2D(pool_size=(2, 2), padding='same', strides=(2,2)))
    model.add(Dropout(0.20))

    ####
    #model.add(Conv2D(128, (3, 3), padding='same', kernel_regularizer=regularizers.l2(0.01),
    #                        bias_initializer='RandomNormal', activation='relu'))
    #model.add(Conv2D(128, (3, 3), padding='same', activation='relu',
    #                        kernel_regularizer=regularizers.l2(0.01),
    #                        bias_initializer='RandomNormal'))
    #model.add(MaxPooling2D(pool_size=(2, 2), padding='same'))
    #model.add(Dropout(0.20))
    ####

    model.add(Flatten())
    model.add(Dense(256, bias_initializer='RandomUniform'))
    model.add(Activation('relu'))
    model.add(Dropout(0.20))
    model.add(Dense(128, bias_initializer='RandomUniform'))
    model.add(Activation('relu'))
    model.add(Dropout(0.20))
    model.add(Dense(1, bias_initializer='RandomUniform'))

    optimizers.Adam(lr=LR, decay=0.0, amsgrad=True)
    model.compile(loss='mean_squared_error', optimizer='adam')

    print("Training model ...")
    history = LossHistory()
    #early_stopping = callbacks.EarlyStopping(monitor='val_loss', patience=10)
    model_info = model.fit(images_train,
                           labels_train,
                           batch_size=BATCH_SIZE,
                           epochs=EPOCHS,
                           validation_data = (images_test, labels_test),
                           verbose=1,
                           callbacks=[history])
                           
    score = model.evaluate(images_test, labels_test, batch_size=16)
    print('SCORE: ', score)

    if SAVE == True:
        # serialize model to JSON
        model_json = model.to_json()
        with open("model_%s.json"%OUTNAME, "w") as json_file:
            json_file.write(model_json)
        # serialize weights to HDF5
        model.save_weights("model_%s.h5"%OUTNAME)
        print("Saved model to disk")
            
            
    #### some plotting
    plt.figure()
    plt.plot(np.arange(len(history.losses)), history.losses)
    plt.xscale('log')
    
    plt.figure(figsize=(5,5))
    plt.scatter(labels_test, model.predict(images_test), alpha=0.5)
    plt.xlabel('True labels')
    plt.ylabel('Predicted labels')
    #plt.plot([-np.pi, np.pi],[-np.pi, np.pi], 'r--')
    
    plt.figure()
    plt.hist(labels_train, bins=100, alpha=0.5, label='Train')
    plt.hist(labels_test, bins=100, alpha=0.5, label='Test')
    plt.hist(model.predict(images_test), bins=100, alpha=0.5, label='Prediction Test')
    plt.title('Labels Distributions')
    plt.legend()
    
    plt.show()
    



if __name__ == '__main__':
    #args = PARSER.parse_args()
    PRRevaluation()
