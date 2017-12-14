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
import matplotlib.patches as patches

import astropy.io.fits as pf
from matplotlib import cm, colors
cmap = cm.get_cmap('viridis')

from PRRutils import build_CNN_tensors



__description__ = 'Data Transformation: from hexagonal to squared pixels to tensors'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-c', '--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('-nevt', '--nevents', type=int, required=False,
                    default=None, help='Number of events to consider')


###################################
####### PATTERN-REC PARAM #########
###################################
PRframe = (38,38)  #Pattern-Recongnition frame


##########################
###### READ THE FILE #####
##########################
def PRRtransform(**kwargs):
    
    f = kwargs['config']
    final_shape = (58,39)
    n_events = kwargs['nevents']
    
    #############################################
    ##### TRANSFORMATION and RETURN TENSORS #####
    #############################################
    f_images = f.replace('.fits', '_images.pkl')
    f_labels = f.replace('.fits', '_labels.pkl')
    if not os.path.exists(f_images):
        images, labels = build_CNN_tensors(f, frame=PRframe, shape=final_shape, nevents=n_events)
        pic.dump(images, open(f_images,'wb'))
        pic.dump(labels, open(f_labels,'wb'))
    else:
        with open(f_images, 'rb') as ff:
            images = pic.load(ff)
        with open(f_labels, 'rb') as fff:
            labels = pic.load(fff)

    #########################
    ###### SOME DRAWING #####
    #########################
    plt.figure()
    plt.imshow(images[1], cmap=cmap)
    plt.title('MC energy = %.2f KeV, MC phi = %.2f rad'%(labels[1][0], labels[1][1]))
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    plt.show()




if __name__ == '__main__':
    args = PARSER.parse_args()
    PRRtransform(**args.__dict__)
