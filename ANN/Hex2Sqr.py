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
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from itertools import product
import astropy.io.fits as pf
from matplotlib import cm, colors
cmap = cm.get_cmap('viridis')


###################################
######### GPD PARAMETERS ##########
###################################

"""See gpd_coordinates.pdf for further details
"""

gpd_dict = {'Ncol' : 300,
            'Nrow' : 352,
            'pitchcol' : 50.0,
            'pitchrow' : 43.3}

Ncols = gpd_dict['Ncol']
Nrows = gpd_dict['Nrow']
pitchcol = gpd_dict['pitchcol']
pitchrow = gpd_dict['pitchrow']
totlenrow = (Ncols - 1/2)*pitchcol
totlencol = (Nrows - 1)*pitchrow

###################################
####### PATTERN-REC PARAM #########
###################################
PRframe = (35,35)  #Pattern-Recongnition frame




def reframe(min_col, max_col, min_row, max_row, frame=PRframe):
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

def buildeventdict(event_params, mc_params):
    """Build the dictionary with all relevant info of an event.
       The reframe hammens here.
        
       ATT: The array of the charges for each pixel, starts fro the bottom-
            left pixel and continues by rows up to the top-right pixel.
            On the other hand the pixel orders starts from the top-lef pixel 
            and continues by columns up to the bottom-right pixel.
    """
    dict = {}
    (en, conv_x, conv_y, pe_en, pe_phi) = mc_params
    (min_col, max_col, min_row, max_row, charge_) = event_params
    totpix = len(charge_)
    charge_ = charge_.reshape((max_row-min_row)+1, (max_col-min_col)+1).T
    charge_ = charge_.reshape(1, totpix)
    max_charge = float(np.amax(charge_))
    colors_ = cmap(charge_[0]/max_charge)
    xpix, ypix = [],[]
    i = np.arange(min_col, max_col+1)
    j = np.arange(min_row, max_row+1)
    ij = list(product(i, j))
    rel_i = np.arange(PRframe[0]+1)
    rel_j = np.arange(PRframe[1]+1)
    rel_ij = list(product(rel_i, rel_j))
    newpix, mcol, Mcol, mrow, Mrow = reframe(min_col, max_col, min_row, max_row)
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
                       colors_[index[0]]]
        else:
            dict[k] = [rel_ij[k], newpix[k], (xpix[k], ypix[k]), cmap(0.)]
    return dict


##########################
###### READ THE FILE #####
##########################
f = '../../sim.fits'
events, mc_energy, mc_abs_x, mc_abs_y, mc_pe_energy, mc_pe_phi = \
                                                          readsimfitsfile(f)
for id, e in enumerate(events[:5]):
    mc_params = (mc_energy[id], mc_abs_x[id], mc_abs_y[id], mc_pe_energy[id], mc_pe_phi[id])
    event_params = (e[5], e[6], e[7], e[8], e[11])
    newpix, min_col, max_col, min_row, max_row = reframe(e[5], e[6], e[7], e[8])
    dict = buildeventdict(event_params, mc_params)
    
    
    #########################
    ###### SOME DRAWING #####
    #########################
    fig = plt.figure(figsize=(12.4, 7.))
    title = 'MC ENERGY = %.2f KeV\nMC X = %.2f mm\nMC Y = %.2f mm'\
                %(mc_energy[id], mc_abs_x[id], mc_abs_y[id])+\
            '\nMC PE ENERGY = %.2f KeV\nMC PE PHI = %.3f'%(mc_pe_energy[id], mc_pe_phi[id])
    grid = plt.GridSpec(4, 6, hspace=0.8, wspace=1.2)
    gpd_ax = fig.add_subplot(grid[:, :-2])
    hex_ax = fig.add_subplot(grid[:-2, -2:])
    sqr_ax = fig.add_subplot(grid[-2:, -2:])
    for k,v in dict.items():
        hexagon = patches.RegularPolygon((v[2][0],v[2][1]), 6,
                                         radius=pitchcol/2000,
                                         color=v[3], orientation=np.pi)
        gpd_ax.add_patch(hexagon)

    gpd_ax.annotate(title,
                    xy= (-6.5, 4.5),
                    xytext=(-6.5, 4.5),
                    size = 12)
    gpd_ax.set_xlim(-totlencol/2000, totlencol/2000)
    gpd_ax.set_ylim(-totlenrow/2000, totlenrow/2000)
    gpd_ax.set_xlabel('[mm]')
    gpd_ax.set_ylabel('[mm]')
    gpd_ax.set_title('GPD Area')
    gpd_ax.grid(alpha=0.3)

    for k,v in dict.items():
        hexagon = patches.RegularPolygon((v[2][0],v[2][1]), 6,
                                     radius=pitchcol/2000,
                                     color=v[3], orientation=np.pi)
        hex_ax.add_patch(hexagon)
    xybottleft = ij2xy(min_col, max_row, Ncols, Nrows, pitchcol, pitchrow)
    xytoprigth = ij2xy(max_col, min_row, Ncols, Nrows, pitchcol, pitchrow)
    hex_ax.plot(mc_abs_x[id], mc_abs_y[id], 'x', color='r')
    hex_ax.plot([xybottleft[0], xytoprigth[0]],
            [np.tan(mc_pe_phi[id])*(xybottleft[0]-mc_abs_x[id])+mc_abs_y[id],
             np.tan(mc_pe_phi[id])*(xytoprigth[0]-mc_abs_x[id])+mc_abs_y[id]],
            'r--', linewidth=0.5)
    hex_ax.set_ylim(xybottleft[1]-2*pitchrow/1000,
                    xytoprigth[1]+2*pitchrow/1000)
    hex_ax.set_xlim(xybottleft[0]-2*pitchcol/1000,
                    xytoprigth[0]+2*pitchcol/1000)
    hex_ax.set_xlabel('[mm]', size=8)
    hex_ax.set_ylabel('[mm]', size=8)
    hex_ax.set_title('%i x %i Area - Hexagonal pixels'%(PRframe[0], PRframe[1]),
                     fontsize=9)

    ###################################################
    ##### FROM HEXAGONAL PIXELS TO SQUARE PIXELS ######
    ###################################################
    dict = hexpix2sqrpix(dict, gpd_dict)

    for k,v in dict.items():
        square = patches.RegularPolygon((v[2][0], v[2][1]), 4,
                                        radius=pitchcol/2000,
                                        color=v[3], orientation=np.pi/4)
        sqr_ax.add_patch(square)
    new_conv_x = hexx2sqrx(mc_abs_x[id], mc_abs_y[id], pitchcol, pitchrow)
    sqr_ax.plot(hexx2sqrx(mc_abs_x[id], mc_abs_y[id], pitchcol, pitchrow), mc_abs_y[id], 'x',
                color='r')
    sqr_ax.plot([hexx2sqrx(xybottleft[0],
                           np.tan(mc_pe_phi[id])*(xybottleft[0]-mc_abs_x[id]) +\
                           mc_abs_y[id], pitchcol, pitchrow),
                 hexx2sqrx(xytoprigth[0],
                           np.tan(mc_pe_phi[id])*(xytoprigth[0]-mc_abs_x[id]) + \
                           mc_abs_y[id], pitchcol, pitchrow)],
            [np.tan(mc_pe_phi[id])*(xybottleft[0]-mc_abs_x[id]) + mc_abs_y[id],
             np.tan(mc_pe_phi[id])*(xytoprigth[0]-mc_abs_x[id]) + mc_abs_y[id]],\
                'r--', linewidth=0.5)
    sqr_ax.set_ylabel('Row ID', size=8)
    sqr_ax.set_xlabel('Column ID', size=8)
    sqr_ax.set_ylim(xybottleft[1], xytoprigth[1])
    sqr_ax.set_xlim(xybottleft[0]-10*pitchcol/1000, xytoprigth[0]+10*pitchcol/1000)
    sqr_ax.set_title('%i x %i Area - Squared pixels'%(PRframe[0], PRframe[1]),
                    fontsize=9)
plt.show()


#if __name__ == "__main__":
#    main()
