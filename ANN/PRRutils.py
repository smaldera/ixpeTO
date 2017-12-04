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


"""Put here all the function and classes need to 
   Pattern-Recognition Recon (PRR)
"""

import imp
import numpy as np
import astropy.io.fits as pf
from itertools import product

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()


###################################
######### GPD PARAMETERS ##########
###################################

"""See gpd_coordinates.pdf for further details
    """
get_var_from_file('gpd_param.py')
gpd_dict = data.GPD_DICT

Ncols = gpd_dict['Ncol']
Nrows = gpd_dict['Nrow']
pitchcol = gpd_dict['pitchcol']
pitchrow = gpd_dict['pitchrow']


def reframe(min_col, max_col, min_row, max_row, frame=(35,35)):
    """ Function to fix the pixel area: it is important to give always
        the same number of inputs (number of pixels) to the NN.
    """
    Nrows = max_row-min_row
    Ncols = max_col-min_col
    drows = frame[0]-Nrows
    dcols = frame[1]-Ncols
    if not np.mod(drows,2):
        r_min_row = min_row - int(drows/2)
        r_max_row = max_row + int(drows/2)
    else:
        r_min_row = min_row - int(drows/2)
        r_max_row = max_row + int(drows/2) + 1
    if not np.mod(dcols,2):
        r_min_col = min_col - int(dcols/2)
        r_max_col = max_col + int(dcols/2)
    else:
        r_min_col = min_col - int(dcols/2)
        r_max_col = max_col + int(dcols/2) + 1
    addedrows = np.arange(r_min_row, r_max_row+1)
    addedcols = np.arange(r_min_col, r_max_col+1)
    addedpix = list(product(addedcols, addedrows))
    return addedpix, r_min_col, r_max_col, r_min_row, r_max_row

def hexx2sqrx(x, y, pitchcol, pitchrow):
    """To get from x coordinate of hexagonal pixeling to x coord
        of squared pixeling. (the y coord remains the same)
    """
    return x-0.5*(pitchcol/pitchrow)*y

def hexpix2sqrpix(event_dict, gpd_dict):
    """To get from hexagonal pixeling to squared pixeling.
    """
    pitchcol = gpd_dict['pitchcol']
    pitchrow = gpd_dict['pitchrow']
    for k,v in event_dict.items():
        new_x = hexx2sqrx(v[2][0], v[2][1], pitchcol, pitchrow)
        event_dict[k] = [v[0], v[1], (new_x, v[2][1]), v[3]]
    return event_dict

def readsimfitsfile(file_path):
    """Function to extract parameters from a Monte Carlo fits file.
    """
    data_f = pf.open(file_path)
    data_f.info()
    events = data_f['EVENTS'].data
    mc_params = data_f['MONTE_CARLO'].data
    mc_energy = mc_params['ENERGY']
    mc_abs_x = mc_params['ABS_X']
    mc_abs_y = mc_params['ABS_Y']
    mc_pe_energy = mc_params['PE_ENE']
    mc_pe_phi = mc_params['PE_PHI']
    gpdcharge = data_f['EVENTS'].data['PIX_PHAS']
    return events, mc_energy, mc_abs_x, mc_abs_y, mc_pe_energy, mc_pe_phi

def ij2xy(i,j, Ncols, Nrows, pitchcol, pitchrow):
    """Function to go from i (col) and j (row) to gpd absolute x and y in mm.
    """
    x = (i - (1/2)*(Ncols - 3/2 + np.mod(j,2)))*pitchcol/1000
    y = (0.5*(Nrows - 1) - j)*pitchrow/1000
    return(x,y)

def buildeventdict(event_params, mc_params, frame=(35,35)):
    """Build the dictionary with all relevant info of an event.
        The reframe hammens here.
        
        ATT: 
        1. The array of the charges for each pixel, starts fro the bottom-
           left pixel and continues by rows up to the top-right pixel.
           On the other hand the pixel orders starts from the top-lef pixel
           and continues by columns up to the bottom-right pixel.
        2. The charge values are set to range from 0 to 1 (needed by TensorFlow)
    """
    dict = {}
    (en, conv_x, conv_y, pe_en, pe_phi) = mc_params
    (min_col, max_col, min_row, max_row, charge_) = event_params
    totpix = len(charge_)
    charge_ = charge_.reshape((max_row-min_row)+1, (max_col-min_col)+1).T
    charge_ = charge_.reshape(1, totpix)
    max_charge = float(np.amax(charge_))
    relcharge_ = charge_[0]/max_charge
    xpix, ypix = [],[]
    i = np.arange(min_col, max_col+1)
    j = np.arange(min_row, max_row+1)
    ij = list(product(i, j))
    rel_i = np.arange(frame[0]+1)
    rel_j = np.arange(frame[1]+1)
    rel_ij = list(product(rel_i, rel_j))
    newpix, mcol, Mcol, mrow, Mrow = reframe(min_col, max_col, min_row, max_row, frame=frame)
    for (ii, jj) in newpix:
        (x, y) = ij2xy(ii, jj, Ncols, Nrows, pitchcol, pitchrow)
        xpix.append(x)
        ypix.append(y)
    keys = np.arange(len(newpix))
    for k in keys:
        index = [ind for ind, (i,j) in enumerate(ij) if i == newpix[k][0] \
                 if j == newpix[k][1]]
        if len(index) != 0:
            dict[k] = [rel_ij[k], newpix[k], (xpix[k], ypix[k]),
                       relcharge_[index[0]]]
        else:
            dict[k] = [rel_ij[k], newpix[k], (xpix[k], ypix[k]), 0.]
    return dict


if __name__ == "__main__":
    
    f = '../../sim.fits'
    events, mc_energy, mc_abs_x, mc_abs_y, mc_pe_energy, mc_pe_phi = \
        readsimfitsfile(f)
