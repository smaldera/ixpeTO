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

#converts mdat and root file in a simple ascii file!!!


import os
import numpy

#import matplotlib.pyplot as plt
import ROOT


from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
from gpdswswig.Io import ixpeInputBinaryFile


#FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data',
#                         'test_fe_500evts.mdat')

#FILE_PATH = '/home/maldera/FERMI/Xipe/rec/data/mdat/xpol_2081.mdat'
FILE_PATH =      '/home/maldera/FERMI/Xipe/rec/ixpeSw/gpdsw/bin/sim.mdat'
ROOT_FILE_PATH = '/home/maldera/FERMI/Xipe/rec/ixpeSw/gpdsw/bin/sim.root'


            

if __name__ == "__main__":



     rootFile = ROOT.TFile(ROOT_FILE_PATH,"open")
     simTree2 = rootFile.Get("ixpe2")
     
     zeroSupThreshold=5
     binary_file = ixpeInputBinaryFile(FILE_PATH)
     clustering = ixpeClustering(zeroSupThreshold)
     #binary_file.buildEventTable()
     num_events=binary_file.numEvents()
     #num_events=100000;
     print " n event  = ",num_events

     miofile = open('test.txt','w')    
     i=0
     #for i in range (0, num_events):
     for entry in simTree2:  

            if i%1000 ==0:
               print i,'...'
        
            photonE=entry.PhoEnergy
            absX=entry.AbsPosX     
            absY=entry.AbsPosY 
            phi=entry.PePhi
            theta=entry.PeTheta
            
            
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
            hit=track.hits()
           
 
            n_hits=track.numHits()
            # print "n_hits  = ",n_hits

            header='event_id= '+str(i)+' nPixels= '+str(n_hits)+' E= '+str(photonE)+' x= '+str(absX)+' y= '+str(absY)+' phi= '+str(phi)+ ' theta= '+str(theta)+'\n'
         
                       
            #print header
            miofile.write(header)
            x=numpy.array([0.]*n_hits)
            y=numpy.array([0.]*n_hits)
            adc=numpy.array([0.]*n_hits)
 
            for jj in range (0,n_hits):
                        
                x[jj]=hit[jj].x
                y[jj]=hit[jj].y
                adc[jj]=hit[jj].pulseHeight
                
                miofile.write(str(hit[jj].x)+"  "+str(hit[jj].y)+"  "+str(hit[jj].pulseHeight)+ '\n')  
            
                
            x=numpy.array([0.]*n_hits)
            y=numpy.array([0.]*n_hits)
            adc=numpy.array([0.]*n_hits)
            i=i+1
            
     miofile.close()
     rootFile.Close()
