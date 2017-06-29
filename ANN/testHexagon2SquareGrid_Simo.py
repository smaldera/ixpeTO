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
from array import array


import matplotlib.pyplot as plt
import ROOT
import math

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
   n_hits=track.numHits()
   print "n_hits  = ",n_hits

   x1=-0.99
   x2=1.01
   y1=-1.01
   y2=1.01
   
   nbinsX=int( (x2-x1)/0.01);
   nbinsY=int((y2-y1)/0.01)
   c=ROOT.TCanvas("c","",0)
   c.Divide(2,2)
   c.GetPad(1).SetLogz()
   c.GetPad(2).SetLogz()
   c.GetPad(3).SetLogz()
    
   
   h2 = ROOT.TH2F("h2","",nbinsX,x1,x2,nbinsY,y1,y2)
   h2.GetXaxis().SetRangeUser(-0.8,1)
   h2.GetYaxis().SetRangeUser(-0.6, 0.6)

   nbinsXhex=int( (x2-x1)/0.029);
   nbinsYhex=int((y2-y1)/0.029)
   
   h2hex=ROOT.TH2Poly("h2poly","",nbinsXhex,x1,x2,nbinsYhex,y1,y2)
   h2hex.Honeycomb(-0.99,-1.01,0.029,nbinsXhex ,nbinsYhex);                    
   h2Off=ROOT.TH2F("h2Off","",40,-20,20,40,-20,20)
   
   
   
   x1=-0.99
   x2=1.01
   
   h2prof=ROOT.TProfile2D("h2prof","",int((x2-x1)/0.052),x1,x2, int((y2-y1)/0.052),y1,y2)
   c.cd(1)
   h2hex.Draw("colz")
   c.cd(2)
   h2prof.Draw("colz")

   
   
   x=numpy.array([0.]*n_hits)
   y=numpy.array([0.]*n_hits)
   adc=numpy.array([0.]*n_hits)
 
   
   hit=track.hits()

   L=0.029
   k=math.sqrt(3.)/2.
   
      
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
  
   esa=[ROOT.TPolyLine]*n_hits  

   
   for i in range (0,n_hits):
          print "x = ",(hit[i].x-0.0125)/0.05,"  x int = ", int(round( (hit[i].x-0.0125)/0.05 +0.1,0 ))   , "y= ",(hit[i].y-0.043/2.)/0.043, " y int =",int(round((hit[i].y-0.043/2.)/0.043,0  ) ) ," adc = ",hit[i].pulseHeight
          h2.Fill(xp[i],yp[i], adc[i])
          h2hex.Fill(xp[i],yp[i], adc[i])
          h2prof.Fill(xp[i],yp[i], adc[i])

          x_int= int(round( (hit[i].x-0.0125)/0.05 +0.1,0 ))
          y_int=int(round((hit[i].y-0.043/2.)/0.043,0  ) )
          h2Off.Fill(x_int,y_int,adc[i])
          
          
          b=array('f',[L+y[i],L/2.+y[i],-L/2.+y[i],-L+y[i],-L/2+y[i], L/2+y[i],L+y[i] ])
          a=array('f',[0+x[i],k*L+x[i],k*L+x[i],0+x[i],-k*L+x[i],-k*L+x[i],0+x[i]]   )
          esa[i]=ROOT.TPolyLine(7,a,b)
          c.cd(1)
          esa[i].Draw("samel")   
          c.cd(2)
          esa[i].Draw("samel")   
          
          
   #c=ROOT.TCanvas("c","",0)
   c.cd(1)      
   #h2.Draw("colZsame")
   #h2prof.Draw("AXIG,same")
   c.cd(2)      
   #h2prof.Draw("colZsame")
  # h2prof.Draw("AXIG,same")
   c.cd(3)
   h2Off.Draw("colZ")
  
   c.Update()
   valore = raw_input('continue?')
   h2.Reset()
   h2prof.Reset()
   #h2hex.Reset()
   

   
            

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

