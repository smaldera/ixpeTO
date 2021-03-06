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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

import astropy.io.fits as pf
from matplotlib import cm, colors
cmap = cm.get_cmap('viridis')

from PRRutils import readsimfitsfile
from PRRutils import buildeventdict
from PRRutils import hexpix2sqrpix
from PRRutils import reframe
from PRRutils import ij2xy
from PRRutils import hexx2sqrx
from PRRutils import get_hex_vertices
from PRRutils import get_sqr_vertices
from PRRutils import complete_square_grid

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
totlenrow = (Ncols - 1/2)*pitchcol
totlencol = (Nrows - 1)*pitchrow

###################################
####### PATTERN-REC PARAM #########
###################################
PRframe = (38, 38)  #Pattern-Recongnition frame


##########################
###### READ THE FILE #####
##########################
f = '../sim.fits'
#f = '../sim500000_line4_pol45.fits'
#f = '../sim_5keV_nopol.fits'
#f = '../sim_Fe55_nopol.fits'
#f = '../sim_line59_nopol.fits'

events, mc_energy, mc_abs_x, mc_abs_y, mc_pe_energy, mc_pe_phi = \
                                                          readsimfitsfile(f)

##########################
##### TRANSFORMATION #####
##########################
for id, e in enumerate(events[:2]):
    mc_params = (mc_energy[id], mc_abs_x[id], mc_abs_y[id], mc_pe_energy[id], mc_pe_phi[id])
    event_params = (e[5], e[6], e[7], e[8], e[11])
    newpix, min_col, max_col, min_row, max_row = reframe(e[5], e[6], e[7], e[8], frame=PRframe)
    dict = buildeventdict(event_params, mc_params, frame=PRframe)
    
    index_ = np.array([item[0] for item in dict.items()])
    pix_centers = np.array([dict[k][2] for k in index_])
    pix_values = np.array([dict[k][3] for k in index_])
    pix_verts = get_hex_vertices(pix_centers, pitchcol/2000.)

    dict = complete_square_grid(dict, frame=PRframe)
    dict = hexpix2sqrpix(dict, gpd_dict)
    index_ = np.array([item[0] for item in dict.items()])
    sqr_centers = np.array([dict[k][2] for k in index_])
    sqr_values = np.array([dict[k][3] for k in index_])
    sqr_verts = get_sqr_vertices(sqr_centers, pitchcol/2000.)
    
    #########################
    ###### SOME DRAWING #####
    #########################
    fig = plt.figure(figsize=(12.4, 7.))
    title = 'MC ENERGY = %.2f KeV\nMC X = %.2f mm\nMC Y = %.2f mm'\
                %(mc_energy[id], mc_abs_x[id], mc_abs_y[id])+\
            '\nMC PE ENERGY = %.2f KeV\nMC PE PHI = %.3f'%(mc_pe_energy[id], mc_pe_phi[id])
    topleftpix = dict[0]
    bottomrightpix = dict[len(dict)-1]
    xybottleft = ij2xy((min_col, max_row), Ncols, Nrows, pitchcol, pitchrow)
    xytoprigth = ij2xy((max_col, min_row), Ncols, Nrows, pitchcol, pitchrow)
    grid = plt.GridSpec(4, 6, hspace=0.8, wspace=1.2)
    gpd_ax = fig.add_subplot(grid[:, :-2])
    hex_ax = fig.add_subplot(grid[:-2, -2:])
    sqr_ax = fig.add_subplot(grid[-2:, -2:])
    hexagon = PolyCollection(pix_verts, array=pix_values, cmap=cmap, edgecolors='none')
    gpd_ax.add_collection(hexagon)
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

    hexagon = PolyCollection(pix_verts, array=pix_values, cmap=cmap, edgecolors='none')
    hex_ax.add_collection(hexagon)
    hex_ax.plot(mc_abs_x[id], mc_abs_y[id], 'x', color='r')
    hex_ax.plot([xybottleft[0], xytoprigth[0]],
            [np.tan(mc_pe_phi[id])*(xybottleft[0]-mc_abs_x[id])+mc_abs_y[id],
             np.tan(mc_pe_phi[id])*(xytoprigth[0]-mc_abs_x[id])+mc_abs_y[id]],
            'r--', linewidth=0.8)
    hex_ax.set_facecolor('black')
    hex_ax.set_ylim(xybottleft[1]-3*pitchrow/1000,
                    xytoprigth[1]+2*pitchrow/1000)
    hex_ax.set_xlim(xybottleft[0]-2*pitchcol/1000,
                    xytoprigth[0]+3*pitchcol/1000)
    hex_ax.set_xlabel('[mm]', size=8)
    hex_ax.set_ylabel('[mm]', size=8)
    hex_ax.set_title('%i x %i Area - Hexagonal pixels'%(PRframe[0], PRframe[1]),
                     fontsize=9)

    ###################################################
    ##### FROM HEXAGONAL PIXELS TO SQUARE PIXELS ######
    ###################################################

    square = PolyCollection(sqr_verts, array=sqr_values, cmap=cmap, edgecolors='none')
    sqr_ax.add_collection(square)
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
                'r--', linewidth=0.8)
    sqr_ax.set_facecolor('black')
    sqr_ax.set_ylabel('Row ID', size=8)
    sqr_ax.set_xlabel('Column ID', size=8)
    sqr_ax.set_ylim(bottomrightpix[2][1]-2*pitchcol/1000, topleftpix[2][1]+2*pitchcol/1000,)
    sqr_ax.set_xlim(topleftpix[2][0]-5*pitchcol/1000, bottomrightpix[2][0]+45*pitchcol/1000)
    sqr_ax.set_title('%i x %i Area - Squared pixels'%(PRframe[0], PRframe[1]),
                    fontsize=9)
plt.show()

