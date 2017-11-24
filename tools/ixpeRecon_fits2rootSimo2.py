# script convertire i .fits in un tree di root (lo so che e' brutto e lento....)
# ---> python  read_fits_hist_root.py  nomefile_recon.fits nomeOutFile.root

import astropy.io.fits as pyfits
import numpy as np
from array import array
import ROOT
import sys


from astropy.io import fits


if __name__ == '__main__':


  
#if number of args not correct
  if (len(sys.argv) != 3):
    print "    usage: python read_fits_hist_root.p  fileName.fits  OutfileName.root"
    sys.exit()

  print "processing file ",sys.argv[1]
  print "out root file ",sys.argv[2]

#get FITS file name
  dataFile = open(sys.argv[1],'r')
  fileName = dataFile.name    
  outFileName=sys.argv[2]



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

  
  f = ROOT.TFile(outFileName, 'recreate')
  
  reconT = ROOT.TTree("reconT", "reconT")
  #bufferID = array('i', [0])
  n = array('i',[0])

  n[0]= len(pha_data['TRK_PHA'])
  #n[0]= 2400
  nEv=n[0] 
  



  """
  TRG_ID= array('i', [0]*n[0])
  TIME = array('d', [0.]*n[0])
  ROI_SIZE = array('i', [0]*n[0])
  MIN_COL = array('i', [0]*n[0])
  MAX_COL = array('i', [0]*n[0])
  MIN_ROW = array('i', [0*n[0]])
  MAX_ROW = array('i', [0]*n[0])
  PAKTNUMB = array('i', [0]*n[0])
  ERR_SUM = array('i', [0]*n[0])
  #PIX_PHAS=?????????????????????????????????
  NUM_PIX = array('i', [0]*n[0])
  TRK_PIX = array('i', [0]*n[0])
  TRK_EFRA = array('d', [0.]*n[0])
  TRK_SN = array('d', [0.]*n[0])
  TRK_BORD =array('i', [0]*n[0])
  
  NUM_CLU = array('i', [0]*n[0])
  TRK_SIZE = array('i', [0]*n[0])
   
  pippoTRK_PHA = array('d', [0.]*nEv)
  
  TRK_PI = array('d', [0.]*n[0])
 
  TRK_PHI1 = array('d', [0.]*n[0])
  TRK_PHI2 = array('d', [0.]*n[0])
 
  TRK_ABSX = array('d', [0.*n[0]])
  TRK_ABSY = array('d', [0.]*n[0])
  TRK_BARX = array('d', [0.]*n[0])
  TRK_BARY = array('d', [0.]*n[0])
  TRK_M2T = array('d', [0.]*n[0])
  TRK_M2L = array('d', [0.]*n[0])
  TRK_M3L = array('d', [0.]*n[0])
  TRK_SKEW = array('d', [0.]*n[0])
  """
  
  TRG_ID= array('i',pha_data['TRG_ID'])
  TIME = array('d',pha_data['TIME'])
  ROI_SIZE = array('i',pha_data['ROI_SIZE'])
  MIN_COL =array( 'i',pha_data['MIN_COL'])
  MAX_COL =array( 'i',pha_data['MAX_COL'])
  MIN_ROW =array('i',pha_data['MIN_ROW'])
  MAX_ROW =array('i', pha_data['MAX_ROW'])
  PAKTNUMB = array('i',pha_data['PAKTNUMB'])
  ERR_SUM = array('i',pha_data['ERR_SUM'])
  #PIX_PHAS=pha_data['PIX_PHAS']?????????????????????????????????
  NUM_PIX =array('i', pha_data['NUM_PIX'])
  #TRK_PIX =array('i', pha_data['TRK_PIX'])
  TRK_EFRA=array('d',pha_data['TRK_EFRA'])
  TRK_SN =array( 'd',pha_data['TRK_SN'])
  TRK_BORD =array('i',pha_data['TRK_BORD'])
  TRK_SIZE =array('i',pha_data['TRK_SIZE'])
  NUM_CLU =array('i',pha_data['NUM_CLU'])
  TRK_PHA =array('d', pha_data['TRK_PHA'])
  TRK_PI = array('d',pha_data['TRK_PI'])
  TRK_PHI2 =array('d',pha_data['TRK_PHI2'])
  TRK_PHI1 =array('d',pha_data['TRK_PHI1'])
  TRK_ABSX =array('d',pha_data['TRK_ABSX'])
  TRK_ABSY =array('d',pha_data['TRK_ABSY'])
  TRK_BARX =array('d',pha_data['TRK_BARX'])
  TRK_BARY =array('d',pha_data['TRK_BARY'])
  TRK_M2T =array('d',pha_data['TRK_M2T'])
  TRK_M2L =array('d',pha_data['TRK_M2L'])
  TRK_M3L =array('d',pha_data['TRK_M2L'])
  TRK_SKEW =array('d',pha_data['TRK_SKEW'])


  reconT.Branch('mynum',n,'mynum/I')
  reconT.Branch('TRG_ID', TRG_ID, 'TRG_ID[mynum]/I')
  reconT.Branch('TIME', TIME, 'TIME[mynum]/D')
  reconT.Branch('ROI_SIZE', ROI_SIZE, 'ROI_SIZE[mynum]/I')
  reconT.Branch('MIN_COL', MIN_COL, 'MIN_COL[mynum]/I')
  reconT.Branch('MAX_COL', MAX_COL, 'MAX_COL[mynum]/I')
  reconT.Branch('MIN_ROW', MIN_ROW, 'MIN_ROW[mynum]/I')
  reconT.Branch('MAX_ROW', MAX_ROW, 'MAX_ROW[mynum]/I')
  reconT.Branch('PAKTNUMB', PAKTNUMB, 'PAKTNUMB[mynum]/I') # che e'?
  reconT.Branch('ERR_SUM',ERR_SUM , 'ERR_SUM[mynum]/I') # che e'? 
  reconT.Branch('NUM_PIX', NUM_PIX , 'NUM_PIX[mynum]/I') # che e'?
  #reconT.Branch('TRK_PIX', TRK_PIX , 'TRK_PIX[mynum]/I') # che e'?
  reconT.Branch('TRK_EFRA', TRK_EFRA , 'TRK_EFRA[mynum]/I') # che e'?
  reconT.Branch('TRK_SN', TRK_SN , 'TRK_SN[mynum]/I') # che e'?
  reconT.Branch('TRK_BORD', TRK_BORD , 'TRK_BORD[mynum]/I') # che e'?
  reconT.Branch('NUM_CLU',NUM_CLU, 'NUM_CLU[mynum]/I')
  reconT.Branch('TRK_SIZE',TRK_SIZE , 'TRK_SIZE[mynum]/I')
  reconT.Branch('TRK_PHA',TRK_PHA, 'TRK_PHA[mynum]/D')
  reconT.Branch('TRK_PI', TRK_PI , 'TRK_PI[mynum]/D')
  reconT.Branch('TRK_PHI1',TRK_PHI1,  'TRK_PHI1[mynum]/D')
  reconT.Branch('TRK_PHI2',TRK_PHI2, 'TRK_PHI2[mynum]/D') 
  reconT.Branch('TRK_ABSX',TRK_ABSX , 'TRK_ABSX[mynum]/D')
  reconT.Branch('TRK_ABSY',TRK_ABSY , 'TRK_ABSY[mynum]/D')
  reconT.Branch('TRK_BARX',TRK_BARX , 'TRK_BARX[mynum]/D')
  reconT.Branch('TRK_BARY',TRK_BARY , 'TRK_BARY[mynum]/D')
  reconT.Branch('TRK_M2T',TRK_M2T  , 'TRK_M2T[mynum]/D')
  reconT.Branch('TRK_M2L',TRK_M2L  , 'TRK_M2L[mynum]/D')
  reconT.Branch('TRK_M3L',TRK_M3L  , 'TRK_M3L[mynum]/D')
  reconT.Branch('TRK_SKEW',TRK_SKEW ,'TRK_SKEW[mynum]/D')
  
  
  reconT.Fill()
  print "TRK_PHA=",TRK_PHA
 
  print "n=",n
    
  reconT.Write()
  f.Close()

  
