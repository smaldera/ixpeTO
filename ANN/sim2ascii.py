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
import sys
import numpy

#import matplotlib.pyplot as plt
import ROOT


from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
from gpdswswig.Io import ixpeInputBinaryFile


#FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data',
#                         'test_fe_500evts.mdat')

#FILE_PATH = '/home/maldera/FERMI/Xipe/rec/data/mdat/xpol_2081.mdat'
FILE_PATH =      '/home/maldera/FERMI/Xipe/rec/ixpeSw/gpdsw/bin/sim_5Kev/sim.mdat'
ROOT_FILE_PATH = '/home/maldera/FERMI/Xipe/rec/ixpeSw/gpdsw/bin/sim_5Kev/sim.root'


            

if __name__ == "__main__":

     if (len(sys.argv) != 4):
          print "    usage: python sim2ascii.py  fileName.mdat filename.root  OutfileName.txt"
          sys.exit()

     print "mdat file= ",sys.argv[1]
     print "root file= ",sys.argv[2]
     print "out file= ",sys.argv[3]

     FILE_PATH=sys.argv[1]
     ROOT_FILE_PATH=sys.argv[2]
     OUT_FILE=sys.argv[3]
     
     rootFile = ROOT.TFile(ROOT_FILE_PATH,"open")
     simTree2 = rootFile.Get("ixpe2")
     
     zeroSupThreshold=5
     binary_file = ixpeInputBinaryFile(FILE_PATH)
     clustering = ixpeClustering(zeroSupThreshold)
     
     num_events=binary_file.numEvents()
     #num_events=100000;
     print " n event  = ",num_events

     hPhi= ROOT.TH1F("hPhi","",80,-4,4) 

     
     miofile = open(OUT_FILE,'w')    
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
            hPhi.Fill(phi)
            
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

            #track recon:
            threshold= zeroSupThreshold
            track.reconstruct(threshold, threshold, False)
            baricenter_X=track.barycenter().x()
            baricenter_Y=track.barycenter().y()
            conversion_point_X=track.absorptionPoint().x()
            conversion_point_Y=track.absorptionPoint().y()
       
            phi0=track.firstPassMomentsAnalysis().phi()
            phi1=track.secondPassMomentsAnalysis().phi()
            pulseH=track.pulseHeight()
 
            n_hits=track.numHits()
            # print "n_hits  = ",n_hits

            header='event_id= '+str(i)+' nPixels= '+str(n_hits)+' E= '+str(photonE)+' x= '+str(absX)+' y= '+str(absY)+' phi= '+str(phi)+ ' theta= '+str(theta)+' x_recStd= '+str(conversion_point_X)+' y_recStd= '+str(conversion_point_Y)+' phi_recStd= '+str(phi1)+' trackPH= '+str(pulseH)+'\n'
         
                       
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

     hPhi.Draw()       
     miofile.close()
     rootFile.Close()
     
