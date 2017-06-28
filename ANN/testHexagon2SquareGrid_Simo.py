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
from gpdswswig.Io import ixpeInputBinaryFile


#FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data',
#                         'test_fe_500evts.mdat')

FILE_PATH = '/home/maldera/FERMI/Xipe/rec/data/mdat/xpol_2081.mdat'



def rec_and_draw(track):
   threshold=5
   #track.reconstruct(threshold, threshold, False) 
   
   hit=track.hits()
   print hit
 
   n_hits=track.numHits()
   print "n_hits  = ",n_hits

   x1=-1
   x2=1
   y1=-1
   y2=1
   
   nbinsX=int( (x2-x1)/0.01);
   nbinsY=int((y2-y1)/0.01)
 
   h2 = ROOT.TH2F("h2","",nbinsX,x1,x2,nbinsY,y1,y2)
   h2.GetXaxis().SetRangeUser(-0.8,1)
   h2.GetYaxis().SetRangeUser(-0.6, 0.6)

   nbinsXhex=int( (x2-x1)/0.05);
   nbinsYhex=int((y2-y1)/0.05)
   
   h2hex=ROOT.TH2Poly("h2poly","",nbinsXhex,x1,x2,nbinsYhex,y1,y2)
   h2hex.Honeycomb(-1,-1,0.05,nbinsXhex ,nbinsYhex);                    

   
   h2prof=ROOT.TProfile2D("h2prof","",nbinsXhex,x1,x2,nbinsYhex,y1,y2)
   


   
   x=numpy.array([0.]*n_hits)
   y=numpy.array([0.]*n_hits)
   adc=numpy.array([0.]*n_hits)
 
   
   hit=track.hits()
        
   for i in range (0,n_hits):
      x[i]=hit[i].x
      y[i]=hit[i].y
      adc[i]=hit[i].pulseHeight

   x0=0
   y0=0
   phi=0*ROOT.TMath.DegToRad()
        
   dx = (x - x0)
   dy = (y - y0)
   xp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
   yp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy    
  
       

   
   for i in range (0,n_hits):
           #print "x = ",hit[i].x, "y= ",hit[i].y," adc = ",hit[i].pulseHeight
         
          h2.Fill(xp[i],yp[i], adc[i])
          #h2hex.Fill(xp[i],yp[i], adc[i])
          #h2prof.Fill(xp[i],yp[i], adc[i])
         
          
   c=ROOT.TCanvas("c","",0)        
   h2.Draw("colZ")
   c.Update()
   valore = raw_input('continue?')
   h2.Reset()
 

   
            

if __name__ == "__main__":

     zeroSupThreshold=5
     binary_file = ixpeInputBinaryFile(FILE_PATH)
     clustering = ixpeClustering(zeroSupThreshold)
     #binary_file.buildEventTable()
     #num_events=binary_file.numEvents()
     num_events=1000;
     print " n event  = ",num_events

        
     for i in range (0, num_events):
            try:    
                evt = binary_file.next()
            except RuntimeError  as  e:
                #print "AAAAAAAAAGGGGHHHHHH!!!!!  e.Value=",str(e)
                if str(e)=='Header mismatch':
                        continue
                else:
                    break         

            tracks = clustering.dbScan(evt)    
            if len(tracks)==0:     # escludo eventi con 0 cluster!!!!
               continue

            track = tracks[0]
            rec_and_draw(track)

