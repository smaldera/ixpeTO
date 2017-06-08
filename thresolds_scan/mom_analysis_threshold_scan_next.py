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

#usage: python -i mom_analysis_threshold_scan.py infile_path [zero-suppression_threshold, default = 5]

#python -i mom_analysis_threshold_scan_classe.py /home/lorenzo/Tesi/data/MDAT_files/001_0001318_data.mdat



import os
import numpy as np
import array
from itertools import product

import matplotlib.pyplot as plt

import ROOT
from ROOT import *

from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
from gpdswswig.Io import ixpeInputBinaryFile


class thresholdScan:

    def __init__(self):

        self.zeroSupThreshold = 5

        self.thresholds1 = array.array('i',(i for i in range(5,21)))    # thresholds are integers (ADC)
        self.thresholds2 = array.array('i',(i for i in range(5,21)))
        self.lenThresholds1 = len(self.thresholds1)
        self.lenThresholds2 = len(self.thresholds2)        


        #Initialize name matrix
        a = np.arange(5, 21)                        # ---> per thresholds1
        b = np.arange(5, 21)                        # ---> per thresholds2
        c = np.array(['%i,%i'%(i,j) for i,j in list(product(self.thresholds1,self.thresholds2))])
        self.name_matrix = np.reshape(c, (16,16))


        #Initialize histograms matrix
        self.hist_array0=16*[ROOT.TH1F]
        self.hist_array1=16*[ROOT.TH1F]
        self.hist_array2=16*[ROOT.TH1F]
        self.hist_array3=16*[ROOT.TH1F]
        self.hist_array4=16*[ROOT.TH1F]
        self.hist_array5=16*[ROOT.TH1F]
        self.hist_array6=16*[ROOT.TH1F]
        self.hist_array7=16*[ROOT.TH1F]
        self.hist_array8=16*[ROOT.TH1F]
        self.hist_array9=16*[ROOT.TH1F]
        self.hist_array10=16*[ROOT.TH1F]
        self.hist_array11=16*[ROOT.TH1F]
        self.hist_array12=16*[ROOT.TH1F]
        self.hist_array13=16*[ROOT.TH1F]
        self.hist_array14=16*[ROOT.TH1F]
        self.hist_array15=16*[ROOT.TH1F]
    
        for i in range(0,16):
            self.hist_array0[i] = ROOT.TH1F(self.name_matrix[0][i],self.name_matrix[0][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array1[i] = ROOT.TH1F(self.name_matrix[1][i],self.name_matrix[1][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array2[i] = ROOT.TH1F(self.name_matrix[2][i],self.name_matrix[2][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array3[i] = ROOT.TH1F(self.name_matrix[3][i],self.name_matrix[3][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array4[i] = ROOT.TH1F(self.name_matrix[4][i],self.name_matrix[4][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array5[i] = ROOT.TH1F(self.name_matrix[5][i],self.name_matrix[5][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array6[i] = ROOT.TH1F(self.name_matrix[6][i],self.name_matrix[6][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array7[i] = ROOT.TH1F(self.name_matrix[7][i],self.name_matrix[7][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array8[i] = ROOT.TH1F(self.name_matrix[8][i],self.name_matrix[8][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array9[i] = ROOT.TH1F(self.name_matrix[9][i],self.name_matrix[9][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array10[i] = ROOT.TH1F(self.name_matrix[10][i],self.name_matrix[10][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array11[i] = ROOT.TH1F(self.name_matrix[11][i],self.name_matrix[11][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array12[i] = ROOT.TH1F(self.name_matrix[12][i],self.name_matrix[12][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array13[i] = ROOT.TH1F(self.name_matrix[13][i],self.name_matrix[13][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array14[i] = ROOT.TH1F(self.name_matrix[14][i],self.name_matrix[14][i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array15[i] = ROOT.TH1F(self.name_matrix[15][i],self.name_matrix[15][i],70,-TMath.Pi(),TMath.Pi()) #80,-4,4

        self.hist_matrix = []

        self.hist_matrix.append(self.hist_array0)
        self.hist_matrix.append(self.hist_array1)
        self.hist_matrix.append(self.hist_array2)
        self.hist_matrix.append(self.hist_array3)
        self.hist_matrix.append(self.hist_array4)
        self.hist_matrix.append(self.hist_array5)
        self.hist_matrix.append(self.hist_array6)
        self.hist_matrix.append(self.hist_array7)
        self.hist_matrix.append(self.hist_array8)
        self.hist_matrix.append(self.hist_array9)
        self.hist_matrix.append(self.hist_array10)
        self.hist_matrix.append(self.hist_array11)
        self.hist_matrix.append(self.hist_array12)
        self.hist_matrix.append(self.hist_array13)
        self.hist_matrix.append(self.hist_array14)
        self.hist_matrix.append(self.hist_array15)
    

    #def getLenThresholds2(self):
    #    return self.lenThresholds2


    #Threshold scan method
    #def threshold_scan(self, zeroSupThreshold=5, num_events=7000):
    #def threshold_scan(self, num_events=7000):
    def threshold_scan(self, zeroSupThreshold = 5, num_events = 7000):

        binary_file = ixpeInputBinaryFile(FILE_PATH)
        #numEventsInFile =  binary_file.numEvents()
        #print numEventsInFile

        print 'zero suppression threshold = ', zeroSupThreshold
        clustering = ixpeClustering(zeroSupThreshold)

        thresholds1 = self.thresholds1
        thresholds2 = self.thresholds2

                                                    # for k in range (0, num_events):
        for k in range (1, 100000000):      # TRY: while il file non e' finito
            if k%1000 == 0:
                print k     
            #if (k > numEventsInFile):
            #    break

            try:
                evt = binary_file.next()    # readEvent(k)
            except RuntimeError  as  e:
                if str(e)=='Header mismatch':
                        continue
                else:
                    break

            tracks = clustering.dbScan(evt)         # applies dbScan ---> tracks is array of clusters
            if len(tracks) != 0:                    # throws cluster with no pixel above 5 ADC
                track = tracks[0]                   # tracks[0] is the principal cluster of the i-th event
                                                    # pi = pulse invariant, pha = ???

                for i in range(0, len(thresholds1)):
                    for j in range(0,len(thresholds2)):
                        track.reconstruct(thresholds1[i], thresholds2[j], False)
                        phi = track.secondPassMomentsAnalysis().phi()
                        # phiDeg = ixpeMath.radToDeg(phi)
                        # self.hist_matrix[i][j].Fill(phiDeg)
                        self.hist_matrix[i][j].Fill(phi)
                        self.hist_matrix[i][j].SetTitle(self.name_matrix[i][j])
        #print self.hist_matrix[15][15].GetMean()
        print 'END'
    
        f.Write()
    
    
    
def test(zeroSupThreshold):
    print zeroSupThreshold
    thrScan = thresholdScan()
    thrScan.threshold_scan(zeroSupThreshold)
    
    
if __name__ == "__main__":


    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-z', '--zero-suppression', type=int, default=5,
                        help = 'zero-suppression threshold')

    args = parser.parse_args()

    #FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data', 'test_fe_500evts.mdat')
    FILE_PATH = args.infile    

    f = TFile("thr_scan_%s_phi.root" %os.path.basename(FILE_PATH).replace('.mdat',''), "recreate") #NEW

    test(args.zero_suppression)

    f.Close()


