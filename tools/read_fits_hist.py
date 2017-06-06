# script per leggere i file .fits
# e disegnare (con matplotlib) alcuni istogrammi

# usage: python -i progName.py fileName.fits phi_distr


import astropy.io.fits as pyfits
import numpy as np
from numpy import array
from array import array
import sys

#import matplotlib SERVE???????????????????????
import matplotlib.pyplot as plt

from astropy.io import fits


if __name__ == '__main__':

#if number of args not correct
  if (len(sys.argv) != 3):
    print "    usage: python -i progName.py fileName.fits phi_dist"
    sys.exit()


#get FITS file name
  dataFile = open(sys.argv[1],'r')
  fileName = dataFile.name

#open FITS file (memmap = True to prevent RAM storage issues)
  pha_list = fits.open(fileName, memmap=True)

#print info on FITS file content
  pha_list.info()

#print columns name, format, unit
  print(pha_list[1].columns)

#load data into a separate variable
  pha_data = pha_list[1].data

  print pha_data[0]	                # row 0
  print pha_data[0][2]	            # row 0, element 2
  print pha_data[0]['track0_phi']   # row 0, element from col 'track0_phi' --> SLOWER PERFORMANCE because an intermediate row object gets created
  print pha_data['track0_phi']	    # column 'track0_phi'
  print pha_data['track0_phi'][0]   # column 'track0_phi', element 0 --> FASTER PERFORMANCE
  print len(pha_data)	            # number of events



#PHI DISTRIBUTION HISTOGRAM
  if sys.argv[2] == 'phi_distr':
    nBinsPhi = 110
    phi_hist = plt.hist(pha_data['track0_phi'], nBinsPhi, color='r')
    plt.title('phi distribution')
    plt.xlabel('phi')
    plt.ylabel('counts')
#    plt.xlim(0,500)
#    plt.ylim(-10,33.2)
    plt.show()

  if sys.argv[2] == 'pulse_height_distr':
    NBINS = 110
    phi_hist = plt.hist(pha_data['track0_pulse_height'], NBINS)
    plt.title('pulse height distribution')
    plt.xlabel('pulse height')
    plt.ylabel('counts')
    plt.show()

  if sys.argv[2] == 'm3_distr':		#m3 is the skweness
    NBINS = 90
    m3_hist = plt.hist(pha_data['track0_skweness'], NBINS)
    plt.title('skweness (m3) distribution')
    plt.xlabel('m3')
    plt.ylabel('counts')
    plt.show()

  if sys.argv[2] == 'm2trans_distr':
    NBINS = 90
    m3_hist = plt.hist(pha_data['track0_mom2trans'], NBINS)
    plt.title('m2trans distribution')
    plt.xlabel('m2trans')
    plt.ylabel('counts')
    plt.show()

  if sys.argv[2] == 'm2long_distr':
    NBINS = 90
    m3_hist = plt.hist(pha_data['track0_mom2long'], NBINS)
    plt.title('m2long distribution')
    plt.xlabel('m2long')
    plt.ylabel('counts')
    plt.show()

  if sys.argv[2] == 'barY:barX':
    plt.hist2d(pha_data['track0_barycenterX'],
               pha_data['track0_barycenterY'],
               bins=100)                      #(x, y, nbins, norm)
    plt.title('barY versus barX')
    plt.xlabel('barX')
    plt.ylabel('barY')
    #plt.colorbar()
    plt.show()


  pha_list.close()


