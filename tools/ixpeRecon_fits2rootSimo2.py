# script convertire i .fits in un tree di root (lo so che e' brutto e lento....)
# ---> python  read_fits_hist_root.py  test_fe_500evts_recon.fits 

import astropy.io.fits as pyfits
import numpy as np
from array import array
import ROOT
import sys


from astropy.io import fits


if __name__ == '__main__':


  
#if number of args not correct
  if (len(sys.argv) != 2):
    print "    usage: python -i progName.py fileName.fits "
    sys.exit()

  print "processing file ",sys.argv[1]


#get FITS file name
  dataFile = open(sys.argv[1],'r')
  fileName = dataFile.name    #????????????????????


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


  print "len pha = ",len(pha_data)		   # number of events
  print "len pha[tke_phi]= ",len( pha_data['TRK_PHI2'])

  f = ROOT.TFile('outSimo2new.root', 'recreate')
  reconT = ROOT.TTree("reconT", "reconT")
 

  bufferID = array('i', [0])
 
   #array(typecode,[initializer])	#--> N.B.: python: i = unsigned integer, d = double
  					# 	       ROOT:   I = unsigned integer, D = double

  eventID = array('i', [0])  	    # maxn* e' da mettere?
  timestamp = array('d', [0.])
  roi_size = array('i', [0])
  roi_min_column = array('i', [0])
  roi_max_column = array('i', [0])
  roi_min_row = array('i', [0])
  roi_max_row = array('i', [0])
  num_tracks = array('i', [0])
  trk_size = array('i', [0])
  trk_pulse_height = array('i', [0])
  trk_phi = array('d', [0.])
  trk_absorptionX = array('d', [0.])
  trk_absorptionY = array('d', [0.])
  trk_barycenterX = array('d', [0.])
  trk_barycenterY = array('d', [0.])
  trk_mom2trans = array('d', [0.])
  trk_mom2long = array('d', [0.])
  trk_skweness = array('d', [0.])

 


  

  reconT.Branch('bufferID', bufferID, 'bufferID/I')
  reconT.Branch('eventID', eventID, 'eventID/I')
  reconT.Branch('roi_size', roi_size, 'roi_size/D')
  reconT.Branch('roi_min_column', roi_min_column, 'roi_min_column/I')
  reconT.Branch('roi_max_column', roi_max_column, 'roi_max_column/I')
  reconT.Branch('roi_min_row', roi_min_row, 'roi_min_row/I')
  reconT.Branch('roi_max_row', roi_max_row, 'roi_max_row/I')
  reconT.Branch('num_tracks', num_tracks, 'num_tracks/I')
  reconT.Branch('trk_size', trk_size, 'trk_size/I')
  reconT.Branch('trk_pulse_height', trk_pulse_height, 'trk_pulse_height/I')
  reconT.Branch('trk_phi', trk_phi, 'trk_phi/D')
  reconT.Branch('trk_absorptionX', trk_absorptionX, 'trk_absorptionX/D')
  reconT.Branch('trk_absorptionY', trk_absorptionY, 'trk_absorptionY/D')
  reconT.Branch('trk_barycenterX', trk_barycenterX, 'trk_barycenterX/D')
  reconT.Branch('trk_barycenterY', trk_barycenterY, 'trk_barycenterY/D')
  reconT.Branch('trk_mom2trans', trk_mom2trans, 'trk_mom2trans/D')
  reconT.Branch('trk_mom2long', trk_mom2long, 'trk_mom2long/D')
  reconT.Branch('trk_skweness', trk_skweness, 'trk_skweness/D')





  
#FILL the TREE
  treeFill=reconT.Fill
  for i in range(0, len(pha_data)):

    eventID[0] = pha_data['eventID'][i]
    num_tracks[0]= pha_data['NUM_CLU'][i]
    trk_size[0]= pha_data['TRK_SIZE'][i]
    trk_pulse_height[0] = pha_data['TRK_PHA'][i]
    trk_phi[0] = pha_data['TRK_PHI2'][i]
    trk_absorptionX[0] =  pha_data['TRK_ABSX'][i]
    trk_absorptionY[0] = pha_data['TRK_ABSY'][i]

    trk_barycenterX[0] = pha_data['TRK_BARX'][i]
    trk_barycenterY[0] = pha_data['TRK_BARY'][i]
    trk_mom2trans[0] = pha_data['TRK_M2T'][i]
    trk_mom2long[0] = pha_data['TRK_M2L'][i]
    trk_skweness[0] = pha_data['TRK_SKEW'][i]

    if i%1000 ==0:
      print "fill event ",i
    #reconT.Fill()
    treeFill()
    
  reconT.Write()
  f.Close()

