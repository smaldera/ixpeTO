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
import numpy

import matplotlib.pyplot as plt
import ROOT


from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
#from gpdswswig.Io import ixpeInputBinaryFile
from gpdswswig.Io import ixpeLvl1FitsFile
from gpdswswig.Event import ixpeEvent
from  gpdswswig.Geometry import *

#FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data',
#                         'test_fe_500evts.mdat')

#FILE_PATH = '/home/maldera/IXPE/data/realData/xpol_2081.mdat'
FILE_PATH = '/home/maldera/IXPE/data/realData/xpol_2081.fits'



def rec_and_draw(track):
   threshold=5
   track.reconstruct(threshold, threshold, False) 

   #loop sulle hit??? 
   #hits=track.ixpeTrack_hits()
   #print "len(hits) = ",len(hits)
   #print "hits",hits
   #pyIt= track.begin
   #print pyIt
   hit=track.hits()
   print (hit)

   print ("---->>>>>>>>>>>> bary = ",track.barycenter().x())
   
   n_hits=track.numHits()
   print ("n_hits  = ",n_hits)

   x1=-2
   x2=2
   y1=-2
   y2=2
   
   nbinsX=int( (x2-x1)/0.01);
   nbinsY=int((y2-y1)/0.01)
 
   h2 = ROOT.TH2F("h2","",nbinsX,x1,x2,nbinsY,y1,y2)
   h2.GetXaxis().SetRangeUser(-0.8,1)
   h2.GetYaxis().SetRangeUser(-0.6, 0.6)

   for i in range (0,n_hits):
           #print "x = ",hit[i].x, "y= ",hit[i].y," adc = ",hit[i].pulseHeight
         
           h2.Fill(hit[i].x,hit[i].y,hit[i].pulseHeight)


   c=ROOT.TCanvas("c","",0)        
   h2.Draw("colZ")
   c.Update()
   valore = input('continue?')
   h2.Reset()
 

def threshold_scan(zeroSupThreshold=5, num_events=20):
        #binary_file = ixpeInputBinaryFile(FILE_PATH)
        binary_file = ixpeLvl1FitsFile(FILE_PATH)
        min_hits=6
        minDensityPoints=4
        # default values= zeroSupThreshold, minTrackHits=6, minDensityPoints=4):
        clustering = ixpeClustering(zeroSupThreshold,min_hits,minDensityPoints)
        #binary_file.buildEventTable()
        #num_events=binary_file.numEvents()
        print (" n event  = ",num_events)

        
        for i in range (0, num_events):
            try:    
                digiEvt = binary_file.next()
            except RuntimeError  as  e:
                #print "AAAAAAAAAGGGGHHHHHH!!!!!  e.Value=",str(e)
                if str(e)=='Header mismatch':
                        continue
                else:
                    break         
            evt = ixpeEvent(digiEvt)

            clustering.dbScan(evt)
            print ('Number of above threshold pixels: %d' %evt.numAboveThresholdPixels())
            print ('Number of noise pixels: %d' % evt.numNoisePixels())
            tracks = ixpeTrackBuilder.buildTracks(evt)
            
           
            if len(tracks)==0:     # escludo eventi con 0 cluster!!!!
               continue


            track = tracks[0]
            rec_and_draw(track)
            threshold=5
            track.reconstruct(threshold, threshold, False)      # la soglia deve essere un intero (ADC)
            
            #phi = track.secondPassMomentsAnalysis().phi()
            #phiDeg = ixpeMath.radToDeg(phi)
            #print "i = ",i," phi = ",phiDeg

            

if __name__ == "__main__":
        threshold_scan(9 ,20)
    
    

