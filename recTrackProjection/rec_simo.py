#!/usr/bin/env python
# Copyright (C) 2007--2016 the X-ray Polarimetry Explorer (XPE) team.
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


import ROOT
import numpy
from array import array
import math
import time
from gpdswpy.matplotlib_ import plt


import smoothing_passabassoSimo  as smooth_simo
from xpeSimo_ttree import *
#from xpeSimo_new import *
from xpeSimo_newMC2 import *


from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
#from gpdswswig.Io import ixpeInputBinaryFile
from gpdswswig.Io import ixpeLvl1FitsFile
from gpdswswig.Event import ixpeEvent
from gpdswswig.MonteCarlo import *
from gpdswpy.event import plot_event
#import gc
from gpdswswig.Calib import ixpeCalibrationSvc, ixpeEventCalibrator


def test(filePath, num_events,raggioCut, dividiBins, baryPadding, findMaxAlg,  zeroSupThreshold=5,  pcubo=0, maxnP=4, Psigma=2, Pthr=0.0001, draw=0):
       """
       """
       print ("file = ",filePath)
       print ("raggioCut = ",raggioCut)
       print ("dividiBins= ",dividiBins)
       print ("baryPadding = ",baryPadding) 
       print ("findMaxAlg = ",findMaxAlg)
       print ("pcubo = ", pcubo)
       print ("maxnP = ", maxnP)
       print ("Psigma ",Psigma)
       print ("Pthr = ",Pthr)
       print ("draw = ",draw)
       print ("zeroSupThresold = ",zeroSupThreshold)
      

       
       h_phi1=ROOT.TH1F("h_phi1","",360,0,360)
       h_phi_tang=ROOT.TH1F("h_phi1_tang","",360,0,360)


       h_x=ROOT.TH1F("h_x","",200,-0.5,0.5)
       h_x1=ROOT.TH1F("h_x1","",200,-0.5,0.5)

       h_y=ROOT.TH1F("h_y","",200,-0.5,0.5)
       h_y1=ROOT.TH1F("h_y1","",200,-0.5,0.5)

       x1=-2 
       x2=4 
       nbinsX=int( (x2-x1)/0.01);
       h1Lcum=ROOT.TH1F("h1L","",int (nbinsX/dividiBins),x1,x2)
       
       #if salva tree?
       event_id=0
       outRootFile=ROOT.TFile("prova1.root","recreate")
       myTree=myTTree()  # from xpeSimo_ttree.py
                     
       #binary_file = ixpeInputBinaryFile(filePath) # questo legge i files mdat ormai obsoleto???
       binary_file= ixpeLvl1FitsFile(filePath)  # questo legge i files fits
       

       # default values= zeroSupThreshold, minTrackHits=6, minDensityPoints=4):
       min_hits=6
       minDensityPoints=4
       clustering = ixpeClustering(zeroSupThreshold,min_hits,minDensityPoints)           
       ixpeCalibrationSvc.setCoherentNoiseOffset(0)
       ixpeCalibrationSvc.setTriggerMiniclusterOffset(0)

       # xpeSimoAA=xpeSimo(raggioCut,dividiBins,baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw)
      # xpeSimoAA.outRootFile= outRootFile
       #print ("event id II= ",xpeSimoAA.event_id)

       for i in range(num_events):
           try:    
                  digiEvt = binary_file.next()
           except StopIteration:
                   #print ("!!!!!!!!!!!!!!!!!! STOP ITERATION EXCEPTION, break loop")
                   break
           except   RuntimeError  as  e:
                #print ("AAAAAAAAAGGGGHHHHHH!!!!!  e.Value=",str(e))
                if str(e)=='Header mismatch':
                        continue
                else:
                    break 
           
           evt = ixpeEvent(digiEvt)
           
           mcInfo=0
           try:
                  mcInfo = binary_file.readMcInfo(i+1)
                  print("mc INFO OK")
           except:      
                 mcInfo=-1 


                 
          # print("adcCounts= ",[adc for adc in evt.adcCounts()], )      

                 
           ixpeEventCalibrator.calibrate(evt)
           clustering.dbScan(evt)
           print ("time = ",evt.timestamp())   
          # plot_event(evt, 25)
          # plt.show()
           
           tracks = ixpeTrackBuilder.buildTracks(evt)
           if len(tracks)==0:     # escludo eventi con 0 cluster!!!!
              continue

           track = tracks[0]
           threshold=zeroSupThreshold
           track.reconstruct(threshold, threshold, False)      # la soglia deve essere un intero (ADC)


           xpeSimoAA=xpeSimo(track,raggioCut,dividiBins,baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw)
           xpeSimoAA.outRootFile= outRootFile

           
           
           xpeSimoAA.event_id=event_id
           xpeSimoAA.McInfo=mcInfo        
           xpeSimoAA.track=track
           xMC=mcInfo.absorbtionPointX
           yMC=mcInfo.absorbtionPointY
           
           recSimo=xpeSimoAA.rec_simo()

           h1Lcum.Add(xpeSimoAA.h1L_smoothSimo,1)
           
           if xpeSimoAA.event_id%100==0: 
                  print ("event id = ",xpeSimoAA.event_id)
       
           if draw:
                   xpeSimoAA.draw_simo()
               
           myTree.Fill(xpeSimoAA)
               
                                           
               
           h_phi1.Fill(xpeSimoAA.phi1*ROOT.TMath.RadToDeg() )
           h_phi_tang.Fill(xpeSimoAA.phiTang*ROOT.TMath.RadToDeg() )
           if (xpeSimoAA.xnew!=-100):
                   h_x.Fill(xpeSimoAA.conversion_point_X-xMC)
                   h_x1.Fill(xpeSimoAA.xnew-xMC)
                   h_y.Fill(xpeSimoAA.conversion_point_Y-yMC)
                   h_y1.Fill(xpeSimoAA.ynew-yMC)


           event_id =event_id+1
           del tracks, mcInfo
           del xpeSimoAA

           
       ####################################3    
       c3=ROOT.TCanvas("c3","",2000,1000)
       c3.Divide(1,2)
              
       h_x1.SetLineColor(2)
       h_y1.SetLineColor(2) 
       c3.cd(1)
       h_x.Draw()
       h_x1.Draw("sames")
       c3.cd(2)
       h_y.Draw()
       h_y1.Draw("sames")


       # scrivi outfile 
       miofile = open('miofile.txt','w')   
       miofile.write(str(raggioCut)+ " "+str(dividiBins)+ "  "+str(baryPadding)+"  "+str(findMaxAlg)+" "+str(pcubo)+ "  "+str(maxnP)+"  "+str(Psigma)+" "+str(Pthr)+"  " +str( h_x1.GetRMS() )+" \n")


       outRootFile.cd()
       h_phi1.Write()
       h_phi_tang.Write()
       h_y.Write()
       h_y1.Write()
       h_x.Write()
       h_x1.Write()
       h1Lcum.Write()
       c3.Write()
       #xpeSimoAA.h1LCumul.Write()
       myTree.treeSimo.Write()
       c3.Update()
       input('continue?')    
       outRootFile.Close()
       del c3
       del h_phi1,    h_phi_tang,  h_y, h_y1, h_x,  h_x1

       
       #valore = input('continue?')    
        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=100,
                        help = 'number of events to be processed')
    parser.add_argument('-z', '--zero-suppression', type=int, default=20,
                        help = 'zero-suppression threshold')
    

    parser.add_argument('-r', '--raggioCut', type=float, default=100,
                        help = 'raggio intorno al fit per accetare i pixel da proiettare')


    parser.add_argument('-divbins', '--dividiBins', type=float, default=1,
                        help = 'fattore per n. di bins histo proiettato')


    parser.add_argument('-baryPadding', '--baryPadding', type=float, default=0.5,
                        help = 'limite distanxa funzione da baricentro ')


    parser.add_argument('-findMaxAlg', '--findMaxAlg', type=int, default=2,
                        help = 'algoritmo ricerca picco e auger... 1 TSpectrum, 2->due gauss - 3-> fit gaus + Expo cutoff  - 4->doppia expCutoff   5-> findPeakSimo ') 
     

    parser.add_argument('-pcubo', '--pcubo', type=int, default=0,
                        help = 'se 0 fissa il parametro di 3 grado a zero -> parabola!  ')
    

    parser.add_argument('-d', '--draw', type=int, default=False,
                        help = 'draw (da aggiungere storage  su file)  (True/ False)  ')
    
    parser.add_argument('-maxnP', '--maxnP', type=int, default=4,
                        help = 'max numner of peaks in TSectrum constructor  ')

    parser.add_argument('-Psigma', '--Psigma', type=int, default=2,
                        help = 'sigma in TSectrum search peak   ')

    parser.add_argument('-Pthr', '--Pthr', type=float, default=0.0001,
                        help = 'threshold in TSectrum search peak   ')
    
    args = parser.parse_args()


    
    start_time = time.time()
    test(args.infile, args.num_events, args.raggioCut,  args.dividiBins, args.baryPadding, args.findMaxAlg, args.zero_suppression, args.pcubo, args.maxnP, args.Psigma, args.Pthr,  args.draw )
    print("--- %s seconds ---" % (time.time() - start_time))
