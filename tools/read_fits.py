# script per leggere i file .fits
# usage: python -i read_fits.py path/to/FITSfile.fits output_file_name.root



import astropy.io.fits as pyfits
import numpy as np
from array import array
import ROOT
from ROOT import *
from math import sqrt

import sys

import matplotlib
import matplotlib.pyplot as plt

from astropy.io import fits


if __name__ == '__main__':

#if number of args not correct
  if (len(sys.argv) != 3): #2
    print "    usage: python -i read_fits.py path/to/FITSfile.fits output_file_name.root"
    sys.exit()

#  print sys.argv[2]

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

#TABLE READING TOOLS
#  print pha_data[0]		   # row 0
#  print pha_data[0][2]		   # row 0, element 2
#  print pha_data[0]['track0_phi']  # row 0, element from col 'track0_phi' --> SLOWER PERFORMANCE because an intermediate row object gets created
#  print pha_data['track0_phi']	   # column 'track0_phi'
#  print pha_data['track0_phi'][0]  # column 'track0_phi', element 0 --> FASTER PERFORMANCE
#  print len(pha_data)		   # number of events


  f = TFile(sys.argv[2], 'recreate')
  reconT = TTree("reconT", "reconT")

  bufferID = array('i', [0])
#            array(typecode,[initializer])	#--> N.B.: python: i = unsigned integer, d = double
  						# 	       ROOT:   I = unsigned integer, D = double
  eventID = array('i', [0])
  timestamp = array('d', [0.])
  roi_size = array('i', [0])
  roi_min_column = array('i', [0])
  roi_max_column = array('i', [0])
  roi_min_row = array('i', [0])
  roi_max_row = array('i', [0])
  num_tracks = array('i', [0])
  track0_size = array('i', [0])
  track0_pulse_height = array('i', [0])
  track0_phi = array('d', [0.])
  track0_absorptionX = array('d', [0.])
  track0_absorptionY = array('d', [0.])
  track0_barycenterX = array('d', [0.])
  track0_barycenterY = array('d', [0.])
  track0_mom2trans = array('d', [0.])
  track0_mom2long = array('d', [0.])
  track0_skweness = array('d', [0.])


#COMMENT WITH ARGUMENTS
  reconT.Branch('bufferID', bufferID, 'bufferID/I')
  reconT.Branch('eventID', eventID, 'eventID/I')
  reconT.Branch('roi_size', roi_size, 'roi_size/D')
  reconT.Branch('roi_min_column', roi_min_column, 'roi_min_column/I')
  reconT.Branch('roi_max_column', roi_max_column, 'roi_max_column/I')
  reconT.Branch('roi_min_row', roi_min_row, 'roi_min_row/I')
  reconT.Branch('roi_max_row', roi_max_row, 'roi_max_row/I')
  reconT.Branch('num_tracks', num_tracks, 'num_tracks/I')
  reconT.Branch('track0_size', track0_size, 'track0_size/I')
  reconT.Branch('track0_pulse_height', track0_pulse_height, 'track0_pulse_height/I')
  reconT.Branch('track0_phi', track0_phi, 'track0_phi/D')
  reconT.Branch('track0_absorptionX', track0_absorptionX, 'track0_absorptionX/D')
  reconT.Branch('track0_absorptionY', track0_absorptionY, 'track0_absorptionY/D')
  reconT.Branch('track0_barycenterX', track0_barycenterX, 'track0_barycenterX/D')
  reconT.Branch('track0_barycenterY', track0_barycenterY, 'track0_barycenterY/D')
  reconT.Branch('track0_mom2trans', track0_mom2trans, 'track0_mom2trans/D')
  reconT.Branch('track0_mom2long', track0_mom2long, 'track0_mom2long/D')
  reconT.Branch('track0_skweness', track0_skweness, 'track0_skweness/D')


#FILL the TREE
  for i in range(0, len(pha_data)):
    eventID[0] = pha_data['eventID'][i]
    track0_pulse_height[0] = pha_data['track0_pulse_height'][i]
    track0_phi[0] = pha_data['track0_phi'][i]
    track0_barycenterX[0] = pha_data['track0_barycenterX'][i]
    track0_barycenterY[0] = pha_data['track0_barycenterY'][i]
    track0_mom2trans[0] = pha_data['track0_mom2trans'][i]
    track0_mom2long[0] = pha_data['track0_mom2long'][i]
    track0_skweness[0] = pha_data['track0_skweness'][i]

    reconT.Fill()
    if i%1000 ==0:
      print '-->', i, 'events read'

#WRITE the TREE on the ROOT file
  reconT.Write()

#  reconT.Print()
#  reconT.Scan('track0_pulse_height')

  pha_list.close()

#  f.Write()
  f.Close()
