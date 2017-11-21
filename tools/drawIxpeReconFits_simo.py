# script convertire i .fits in un tree di root (lo so che e' brutto e lento....)
# ---> python  read_fits_hist_root.py  nomefile_recon.fits nomeOutFile.root

import astropy.io.fits as pyfits
import numpy as np
import numpy.ma as ma 

from array import array
import ROOT
import sys


from astropy.io import fits


if __name__ == '__main__':

  print "processing file ",sys.argv[1]
  
#if number of args not correct
  if (len(sys.argv) != 3):
    print "    usage: python read_fits_hist_root.p  fileName.fits  OutfileName.root"
    sys.exit()


  hSize=ROOT.TH1F("h_size","tkr_size",1000,0,1000)
  h_numTraks=ROOT.TH1F("h_numTracks","numTracks",100,0,100)
  h_pulseHeight=ROOT.TH1F("h_pulseHeight","pulseHeight ",1000,0,10000)
  h_pulseInvariant=ROOT.TH1F("h_pulseInvariant","pulseInvariant ",1000,0,10000)
 
  h_phi1=ROOT.TH1F("h_phi1","phi1 ",360,-3.1415,3.1415)  
  h_phi2=ROOT.TH1F("h_phi2","phi2 ",360,-3.1415,3.1415)
  h_absX=ROOT.TH1F("h_absX","absX ",200, -10,10)
  h_absY=ROOT.TH1F("h_absY","absY ",200, -10,10)
  h_baryX=ROOT.TH1F("h_baryX","baryX ",200, -10,10)
  h_baryY=ROOT.TH1F("h_baryY","baryY ",200, -10,10)
  h_absXY=ROOT.TH2F("h_absXY","absX absY ",200, -10,10,200,-10,10)
  h_pulse_size=ROOT.TH2F("h_pulse-size","pulse height vs size ",1000,0,10000,100,0,1000)
 

    
  print "filelist ",sys.argv[1]
  print "out root file ",sys.argv[2]


  myfile=open(sys.argv[1],"r")

  for reconFitsFile in myfile:
    print "processing ",reconFitsFile
    #get FITS file name
    dataFile = open(reconFitsFile[0:len(reconFitsFile)-1],'r')
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
  

    # eventID = array('i', [0])  
    # timestamp = array('d', [0.])
    # roi_size = array('i', [0])
    # roi_min_column = array('i', [0])
    # roi_max_column = array('i', [0])
    # roi_min_row = array('i', [0])
    # roi_max_row = array('i', [0])

    trk_sizeNP=np.array(pha_data['TRK_SIZE'])
    num_tracksNP=np.array( pha_data['NUM_CLU'])
    trk_pulse_heightNP =np.array( pha_data['TRK_PHA'])
    trk_pulse_invNP =np.array( pha_data['TRK_PI'])
    trk_phi2NP =np.array( pha_data['TRK_PHI2'])
    trk_phi1NP =np.array( pha_data['TRK_PHI1'])
    trk_absorptionXNP =np.array(  pha_data['TRK_ABSX'])
    trk_absorptionYNP =np.array( pha_data['TRK_ABSY'])
    trk_barycenterXNP =np.array( pha_data['TRK_BARX'])
    trk_barycenterYNP =np.array( pha_data['TRK_BARY'])
    trk_mom2transNP =np.array(pha_data['TRK_M2T'])
    trk_mom2longNP =np.array( pha_data['TRK_M2L'])
    trk_skwenessNP =np.array( pha_data['TRK_SKEW'])

    ############## CUT #######################3
    
    print   "trk_sizeNP", trk_sizeNP
    trk_sizeNP_masked=ma.MaskedArray(trk_sizeNP, ~(trk_pulse_heightNP<100000 )   )  # la maskera chiede cosa voglio tagliare!!!! per questo la nego!! 
    trk_sizeNP_masked2=ma.MaskedArray(trk_sizeNP, ~(trk_pulse_heightNP>0  )   )  # la maskera chiede cosa voglio tagliare!!!! per questo la nego!! 

    
    
    print " trk_sizeNP_masked = ",trk_sizeNP_masked
    mask1=trk_sizeNP_masked.mask
    mask2=trk_sizeNP_masked2.mask

    #devo fare l'or bit a bit di tutte le maschere!!!!
    mask=mask1 | mask2
    
    print "mask = ",mask
    num_tracksNP_masked=ma.masked_array( num_tracksNP ,mask)
    trk_pulse_heightNP_masked =ma.masked_array( trk_pulse_heightNP,mask)
    trk_pulse_invNP_masked =ma.masked_array(  trk_pulse_invNP,mask)
    trk_phi2NP_masked =ma.masked_array( trk_phi2NP,mask)
    trk_phi1NP_masked=ma.masked_array(  trk_phi1NP,mask)
    trk_absorptionXNP_masked =ma.masked_array( trk_absorptionXNP,mask)
    trk_absorptionYNP_masked =ma.masked_array(  trk_absorptionYNP,mask)
    trk_barycenterXNP_masked=ma.masked_array(   trk_barycenterXNP,mask)
    trk_barycenterYNP_masked=ma.masked_array(  trk_barycenterYNP,mask)
    trk_mom2transNP_masked=ma.masked_array( trk_mom2transNP,mask)
    trk_mom2longNP_masked=ma.masked_array(  trk_mom2longNP,mask)
    trk_skwenessNP_masked =ma.masked_array(  trk_skwenessNP,mask)
    
    print trk_sizeNP_masked[~mask].data


    
    
    trk_size=array('d', trk_sizeNP_masked[~mask].data)
    num_tracks=array('d', num_tracksNP_masked[~mask].data         )
    trk_pulse_height =array('d',  trk_pulse_heightNP_masked[~mask].data )
    trk_pulse_inv =array('d',  trk_pulse_invNP_masked[~mask].data )

    trk_phi2 =array('d',   trk_phi2NP_masked[~mask].data          )
    trk_phi1 =array('d',    trk_phi1NP_masked[~mask].data       )
    trk_absorptionX =array('d',trk_absorptionXNP_masked[~mask].data  )
    trk_absorptionY =array('d', trk_absorptionYNP_masked[~mask].data )
    trk_barycenterX =array('d',   trk_barycenterXNP_masked[~mask].data )
    trk_barycenterY =array('d',   trk_barycenterYNP_masked[~mask].data )
    trk_mom2trans =array('d',    trk_mom2transNP_masked[~mask].data )
    trk_mom2long =array('d',    trk_mom2longNP_masked[~mask].data  )
    trk_skweness =array('d',   trk_skwenessNP_masked[~mask].data    )
    

    w= array('d', [1.]*len(trk_size))
    n=len(trk_size) 
    #histogram of trk_size:
    hSize.FillN(n,trk_size,w)
    #hSize.Draw("hist")
  
    #histogram of num_traks:
    h_numTraks.FillN(n,num_tracks,w)

    #histogram of pulse_Height:
    h_pulseHeight.FillN(n,trk_pulse_height,w)
    
    #histogram of pulse_inv:
    h_pulseInvariant.FillN(n,trk_pulse_inv,w)

    #histogram of phi1:
    h_phi1.FillN(n,trk_phi1,w)
  
    #histogram of phi1:
    h_phi2.FillN(n,trk_phi2,w)
    
    #histogram of absX:
    h_absX.FillN(n,trk_absorptionX,w)
    
    #histogram of absY:
    h_absY.FillN(n,trk_absorptionY,w)
    
    #histogram of baryX:
    h_baryX.FillN(n,trk_barycenterX,w)
    
    #histogram of baryY:
    h_baryY.FillN(n,trk_barycenterY,w)
    
    #x-y rec
    h_absXY.FillN(n,trk_absorptionX,trk_absorptionY ,w)
    #Pheight-size rec
    h_pulse_size.FillN(n,trk_pulse_height,trk_size,w)



  f = ROOT.TFile(outFileName, 'recreate')
    

  hSize.Write()
  h_numTraks.Write()
  h_pulseHeight.Write()
  h_pulseInvariant.Write()

  h_phi1.Write()
  h_phi2.Write()
  h_absX.Write()
  h_absY.Write()
  h_baryX.Write()
  h_baryY.Write()
  h_absXY.Write()
  h_pulse_size.Write()

  f.Close()

