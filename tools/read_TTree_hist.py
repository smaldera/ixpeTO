# script per leggere i TTree nei file ROOT
# e disegnare (con ROOT) istogrammi vari
# usage: python -i read_TTree_hist.py file_name.root


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
  if (len(sys.argv) != 2):
    print "    usage: python -i read_TTree_hist.py file_name.root"
    sys.exit()

# GET input file
  inFile = TFile(sys.argv[1])

#GET TTree from input file
  reconT = inFile.Get('reconT')


  c = TCanvas("c","ixpeRecon",0)
  c.Divide(2,2)


  c.cd(1)
# PHI distribution HISTOGRAM
#  reconT.Draw('track0_phi')
  nBinsPhi = 90
  minPhi = reconT.GetMinimum('track0_phi')
  maxPhi = reconT.GetMaximum('track0_phi')
  hPhi = TH1F('hPhi','distrPhi',nBinsPhi,minPhi-50,maxPhi+50)
  hPhi.GetXaxis().SetTitle('phi (rad)')
  hPhi.GetYaxis().SetTitle('counts')
  reconT.Draw('track0_phi>>hPhi')
  fitF = TF1("fitF", "[0]+[1]*(sin(x+[2])*sin(x+[2]))", minPhi, maxPhi) #fitta la distribuzione
  hPhi.Fit("fitF","MR")
  hPhi.SetLineColor(kRed)


  c.cd(2)
#SKWENESS (m3) DISTRIBUTION histogram
#  reconT.Draw('track0_skweness')
  nBinsM3 = 50
  minM3 = reconT.GetMinimum('track0_skweness')
  maxM3 = reconT.GetMaximum('track0_skweness')
  hM3 = TH1F('hM3','distrM3',nBinsM3,minM3-0.4,maxM3+0.4)
  hM3.GetXaxis().SetTitle('skewness (mm^3)') #??? u.d.m. giusta???
  hM3.GetYaxis().SetTitle('counts')
  reconT.Draw('track0_skweness>>hM3')
  hM3.SetLineColor(kOrange)


  c.cd(3)
#barX DISTRIBUTION histogram
#  reconT.Draw('track0_barycenterX')
  nBinsBarX = 90
  minBarX = reconT.GetMinimum('track0_barycenterX')
  maxBarX = reconT.GetMaximum('track0_barycenterX')
  #minBarX = reconT.GetMinimum('track0_absorptionX')
  #maxBarX = reconT.GetMaximum('track0_absorptionX')
  hBarX = TH1F('hBarX','distrBarX',nBinsBarX,minBarX-0.2,maxBarX+0.2)
  hBarX.GetXaxis().SetTitle('BarX (mm)') #??? u.d.m. giusta???
  hBarX.GetYaxis().SetTitle('counts')
  #reconT.Draw('track0_barycenterX>>hBarX')
  reconT.Draw('track0_barycenterX>>hBarX')
  hBarX.SetLineColor(kAzure)
#  print ("BarX distr: RMS = ", hBarX.GetRMS(), " +- ", hBarX.GetRMSError())


  c.cd(4)
#barY DISTRIBUTION histogram
#  reconT.Draw('track0_barycenterY')
  nBinsBarY = 90
  minBarY = reconT.GetMinimum('track0_barycenterY')
  maxBarY = reconT.GetMaximum('track0_barycenterY')
  hBarY = TH1F('hBarY','distrBarY',nBinsBarY,minBarY-0.2,maxBarY+0.2)
  hBarY.GetXaxis().SetTitle('BarY (mm)') #??? u.d.m. giusta???
  hBarY.GetYaxis().SetTitle('counts')
  reconT.Draw('track0_barycenterY>>hBarY')
  hBarY.SetLineColor(kGreen+2)
  print ("BarY distr: RMS = ", hBarY.GetRMS(), " +- ", hBarY.GetRMSError())

  inFile.Close()




  '''
  OPPURE, simpler way:

  c = TCanvas("c","ixpeRecon",0)
  c.Divide(2,2)

  c.cd(1)
  reconT.Draw('track0_absorptionX')
  c.cd(2)
  reconT.Draw('track0_absorptionY')
  c.cd(3)
  reconT.Draw('track0_pulse_height')
  c.cd(4)
  reconT.Draw('track0_phi')

  '''
