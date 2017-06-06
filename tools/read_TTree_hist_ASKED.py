# script per leggere i TTree nei file ROOT
# e disegnare (con ROOT) l'istogramma richiesto
# usage: python -i read_TTree_hist.py file_name.root requested_output


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
  if (len(sys.argv) != 3):
    print "    usage: python -i read_TTree_hist_ASKED.py file_name.root requested_output"
    sys.exit()

# GET input file
  inFile = TFile(sys.argv[1])

#GET TTree from input file
  reconT = inFile.Get('reconT')

#????
  if sys.argv[2] == 'pulse_height_distr':
    nBinsPulseHeight = 100
    cPulseHeight = TCanvas('cPulseHeight','distrPulseHeight', 200, 10, 600, 400)
    cPulseHeight.cd()
    hPulseHeight = TH1F('hPulseHeight','distrPulseHeight',nBinsPulseHeight,reconT.GetMinimum('track0_pulse_height'),reconT.GetMaximum('track0_pulse_height'))
    hPulseHeight.GetXaxis().SetTitle('pulse height (ADC counts)')
    hPulseHeight.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_pulse_height>>hPulseHeight')
    hPulseHeight.SetLineColor(kBlue)

  '''
# PHI_DISTRIBUTION histogram
  if sys.argv[2] == 'phi_distr':
    nBinsPhi = 80
    #minPhi = reconT.GetMinimum('track0_phi')
    #maxPhi = reconT.GetMaximum('track0_phi')
    cPhi = TCanvas('cPhi','distrPhi', 200, 10, 600, 400)
    cPhi.cd()
#   reconT.Draw("track0_phi")
#    hPhi = TH1F('hPhi','distrPhi',nBinsPhi,reconT.GetMinimum('track0_phi'),reconT.GetMaximum('track0_phi'))
#    hPhi = TH1F('hPhi','distrPhi',nBinsPhi,-200,200)
#    hPhi.GetXaxis().SetTitle('phi (rad)')
#    hPhi.GetYaxis().SetTitle('counts')
#    reconT.Draw('track0_phi>>hPhi')
#    reconT.Draw("track0_phi>>hPhi","","same")
#    hPhi.SetLineColor(kRed)
    hBarX = TH1F('hBarX','distrBarX',80,-200,200)
    hBarX.GetXaxis().SetTitle('BarX (mm)') #??? u.d.m. giusta???
    hBarX.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_phi>>hBarX')
  '''

#phi DISTRIBUTION histogram
  if sys.argv[2] == 'phi_distr':
    nBinsBarX = 90
    minBarX = reconT.GetMinimum('track0_phi')
    maxBarX = reconT.GetMaximum('track0_phi')
    cBarX = TCanvas('cBarX','distrBarX', 200, 10, 600, 400)
    cBarX.cd()
#    hBarX = TH1F('hBarX','distrBarX',nBinsBarX,reconT.GetMinimum('track0_barycenterX'),reconT.GetMaximum('track0_barycenterX'))
    hBarX = TH1F('hBarX','distrBarX',nBinsBarX,minBarX-50,maxBarX+50)
    hBarX.GetXaxis().SetTitle('BarX (mm)') #??? u.d.m. giusta???
    hBarX.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_phi>>hBarX')
#    reconT.Draw('track0_barycenterX>>hBarX','track0_barycenterX>0')
    hBarX.SetLineColor(kRed)
    print ("phi distr: RMS = ", hBarX.GetRMS(), " +- ", hBarX.GetRMSError())

  '''
#SKWENESS (m3) DISTRIBUTION histogram
  if sys.argv[2] == 'm3_distr':		#m3 is skweness
    nBinsM3 = 90
    cM3 = TCanvas('cM3','distrM3', 200, 10, 600, 400)
    cM3.cd()
#    hM3 = TH1F('hM3','distrM3',nBinsM3,reconT.GetMinimum('track0_skweness'),reconT.GetMaximum('track0_skweness'))
    hM3 = TH1F('hM3','distrM3',90,-1,1)
    hM3.GetXaxis().SetTitle('skweness (mm^3)') #??? u.d.m. giusta???
    hM3.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_skweness>>hM3')
    hM3.SetLineColor(kGreen)
  '''


#SKEWNESS DISTRIBUTION histogram
  if sys.argv[2] == 'm3_distr':
    nBinsBarY = 90
    cBarY = TCanvas('cBarY','distrBarY', 200, 10, 600, 400)
    cBarY.cd()
    hBarY = TH1F('hBarY','distrBarY',nBinsBarY,reconT.GetMinimum('track0_skweness'),reconT.GetMaximum('track0_skweness'))
    hBarY.GetXaxis().SetTitle('BarY (mm)') #??? u.d.m. giusta???
    hBarY.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_skweness>>hBarY')
    hBarY.SetLineColor(kRed)
    print ("skewness distr: RMS = ", hBarY.GetRMS(), " +- ", hBarY.GetRMSError())


  '''
#m2trans DISTRIBUTION histogram
  if sys.argv[2] == 'm2trans_distr':
    nBinsM2trans = 90
    cM2trans = TCanvas('cM2trans','distrM2trans', 200, 10, 600, 400)
    cM2trans.cd()
    hM2trans = TH1F('hM2trans','distrM2trans',nBinsM2trans,reconT.GetMinimum('track0_mom2trans'),reconT.GetMaximum('track0_mom2trans'))
    hM2trans.GetXaxis().SetTitle('m2_trans (mm^2)') #??? u.d.m. giusta???
    hM2trans.GetYaxis().SetTitle('counts')
#    reconT.Draw('track0_mom2trans>>hM2trans')
    reconT.Draw('track0_mom2trans>>hM2trans')
    hM2trans.SetLineColor(kViolet)
  '''

#m2trans DISTRIBUTION histogram
  if sys.argv[2] == 'm2trans_distr':
    nBinsBarY = 90
    cBarY = TCanvas('cBarY','distrBarY', 200, 10, 600, 400)
    cBarY.cd()
    hBarY = TH1F('hBarY','distrBarY',nBinsBarY,reconT.GetMinimum('track0_mom2trans')-0.5,reconT.GetMaximum('track0_mom2trans')+0.5)
    hBarY.GetXaxis().SetTitle('BarY (mm)') #??? u.d.m. giusta???
    hBarY.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_mom2trans>>hBarY')
    hBarY.SetLineColor(kRed)
    print ("mom2trans distr: RMS = ", hBarY.GetRMS(), " +- ", hBarY.GetRMSError())


  '''
#m2long DISTRIBUTION histogram
  if sys.argv[2] == 'm2long_distr':
    nBinsM2long = 90
    cM2long = TCanvas('cM2long','distrM2long', 200, 10, 600, 400)
    cM2long.cd()
    hM2long = TH1F('hM2long','distrM2long',nBinsM2long,reconT.GetMinimum('track0_mom2long'),reconT.GetMaximum('track0_mom2long'))
    hM2long.GetXaxis().SetTitle('m2_long (mm^2)') #??? u.d.m. giusta???
    hM2long.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_mom2long>>hM2long')
    hM2long.SetLineColor(kOrange)
  '''


#m2long DISTRIBUTION histogram
  if sys.argv[2] == 'm2long_distr':
    nBinsBarY = 90
    cBarY = TCanvas('cBarY','distrBarY', 200, 10, 600, 400)
    cBarY.cd()
    hBarY = TH1F('hBarY','distrBarY',nBinsBarY,reconT.GetMinimum('track0_mom2long')-0.5,reconT.GetMaximum('track0_mom2long')+0.5)
    hBarY.GetXaxis().SetTitle('BarY (mm)') #??? u.d.m. giusta???
    hBarY.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_mom2long>>hBarY')
    hBarY.SetLineColor(kRed)
    print ("mom2long distr: RMS = ", hBarY.GetRMS(), " +- ", hBarY.GetRMSError())



#barX DISTRIBUTION histogram
  if sys.argv[2] == 'barX_distr':
    nBinsBarX = 90
    minBarX = reconT.GetMinimum('track0_barycenterX')
    maxBarX = reconT.GetMaximum('track0_barycenterX')
    cBarX = TCanvas('cBarX','distrBarX', 200, 10, 600, 400)
    cBarX.cd()
#    hBarX = TH1F('hBarX','distrBarX',nBinsBarX,reconT.GetMinimum('track0_barycenterX'),reconT.GetMaximum('track0_barycenterX'))
    hBarX = TH1F('hBarX','distrBarX',nBinsBarX,minBarX,maxBarX)
    hBarX.GetXaxis().SetTitle('BarX (mm)') #??? u.d.m. giusta???
    hBarX.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_barycenterX>>hBarX')
#    reconT.Draw('track0_barycenterX>>hBarX','track0_barycenterX>0')
    hBarX.SetLineColor(kAzure)
    print ("BarX distr: RMS = ", hBarX.GetRMS(), " +- ", hBarX.GetRMSError())


#barY DISTRIBUTION histogram
  if sys.argv[2] == 'barY_distr':
    nBinsBarY = 90
    cBarY = TCanvas('cBarY','distrBarY', 200, 10, 600, 400)
    cBarY.cd()
    hBarY = TH1F('hBarY','distrBarY',nBinsBarY,reconT.GetMinimum('track0_barycenterY'),reconT.GetMaximum('track0_barycenterY'))
    hBarY.GetXaxis().SetTitle('BarY (mm)') #??? u.d.m. giusta???
    hBarY.GetYaxis().SetTitle('counts')
    reconT.Draw('track0_barycenterY>>hBarY')
    hBarY.SetLineColor(kRed)
    print ("BarY distr: RMS = ", hBarY.GetRMS(), " +- ", hBarY.GetRMSError())


#barX:barY scatter plot
  if sys.argv[2] == 'barY:barX':
    nBinsBarX = 90
    nBinsBarY = 90
    cBarXBarY = TCanvas('cBarXBarY','BarX:BarY', 200, 10, 600, 400)
    cBarXBarY.cd()
    hBarXBarY = TH2F('hBarXBarY','barX:barY',nBinsBarX,reconT.GetMinimum('track0_barycenterX'),reconT.GetMaximum('track0_barycenterX'),nBinsBarY,reconT.GetMinimum('track0_barycenterY'),reconT.GetMaximum('track0_barycenterY'))
    hBarXBarY.GetXaxis().SetTitle('BarX (mm)')
    hBarXBarY.GetYaxis().SetTitle('BarY (mm)')
#    reconT.Draw('track0_barycenterX:track0_barycenterY>>hBarXBarY',"","*")
    reconT.Draw('track0_barycenterY:track0_barycenterX',"","*")   #Draw('y:x')


  inFile.Close()

