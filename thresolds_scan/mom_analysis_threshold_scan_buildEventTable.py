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



#DA USARE QUANDO SI SARÀ RISOLTO IL PROBLEMA DELLA CLASSE ixpesw/Io/ixpeEventBinaryFile

import os
import numpy
import array

import matplotlib.pyplot as plt

import ROOT
from ROOT import *

from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
from gpdswswig.Io import ixpeInputBinaryFile




#FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data', 'test_fe_500evts.mdat')


class thresholdScan:

    def __init__(self):

        self.zeroSupThreshold = 5

        self.thresholds1 = array.array('i',(i for i in range(5,21)))
        self.lenThresholds1 = len(self.thresholds1)
        self.thresholds2 = array.array('i',(i for i in range(5,21)))
        self.lenThresholds2 = len(self.thresholds2)        

        '''
        import numpy as np
        from itertools import product
        a = np.arange(5, 21)    ---> per thresholds1
        >>> #a = np.array([1,2,3])
        >>> b = np.array([1,2,3])
        >>> c = np.array(['%i,%i'%(i,j) for i,j in list(product(a,b))])
        >>> cc = np.reshape(c, (3,3)) ---> la mia name_matrix
        cc[1][1]
        '''

        self.name_h0 = ["1,1","1,2","1,3","1,4","1,5","1,6","1,7","1,8","1,9","1,10","1,11","1,12","1,13","1,14","1,15","1,16"]
        self.name_h1 = ["2,1","2,2","2,3","2,4","2,5","2,6","2,7","2,8","2,9","2,10","2,11","2,12","2,13","2,14","2,15","2,16"]
        self.name_h2 = ["3,1","3,2","3,3","3,4","3,5","3,6","3,7","3,8","3,9","3,10","3,11","3,12","3,13","3,14","3,15","3,16"]
        self.name_h3 = ["4,1","4,2","4,3","4,4","4,5","4,6","4,7","4,8","4,9","4,10","4,11","4,12","4,13","4,14","4,15","4,16"]
        self.name_h4 = ["5,1","5,2","5,3","5,4","5,5","5,6","5,7","5,8","5,9","5,10","5,11","5,12","5,13","5,14","5,15","5,16"]
        self.name_h5 = ["6,1","6,2","6,3","6,4","6,5","6,6","6,7","6,8","6,9","6,10","6,11","6,12","6,13","6,14","6,15","6,16"]
        self.name_h6 = ["7,1","7,2","7,3","7,4","7,5","7,6","7,7","7,8","7,9","7,10","7,11","7,12","7,13","7,14","7,15","7,16"]
        self.name_h7 = ["8,1","8,2","8,3","8,4","8,5","8,6","8,7","8,8","8,9","8,10","8,11","8,12","8,13","8,14","8,15","8,16"]
        self.name_h8 = ["9,1","9,2","9,3","9,4","9,5","9,6","9,7","9,8","9,9","9,10","9,11","9,12","9,13","9,14","9,15","9,16"]
        self.name_h9 = ["10,1","10,2","10,3","10,4","10,5","10,6","10,7","10,8","10,9","10,10","10,11","10,12","10,13","10,14","10,15","10,16"]
        self.name_h10= ["11,1","11,2","11,3","11,4","11,5","11,6","11,7","11,8","11,9","11,10","11,11","11,12","11,13","11,14","11,15","11,16"]
        self.name_h11= ["12,1","12,2","12,3","12,4","12,5","12,6","12,7","12,8","12,9","12,10","12,11","12,12","12,13","12,14","12,15","12,16"]
        self.name_h12= ["13,1","13,2","13,3","13,4","13,5","13,6","13,7","13,8","13,9","13,10","13,11","13,12","13,13","13,14","13,15","13,16"]
        self.name_h13= ["14,1","14,2","14,3","14,4","14,5","14,6","14,7","14,8","14,9","14,10","14,11","14,12","14,13","14,14","14,15","14,16"]
        self.name_h14= ["15,1","15,2","15,3","15,4","15,5","15,6","15,7","15,8","15,9","15,10","15,11","15,12","15,13","15,14","15,15","15,16"]
        self.name_h15= ["16,1","16,2","16,3","16,4","16,5","16,6","16,7","16,8","16,9","16,10","16,11","16,12","16,13","16,14","16,15","16,16"]
    
        self.name_matrix = []
        self.name_matrix.append(self.name_h0)
        self.name_matrix.append(self.name_h1)
        self.name_matrix.append(self.name_h2)
        self.name_matrix.append(self.name_h3)
        self.name_matrix.append(self.name_h4)
        self.name_matrix.append(self.name_h5)
        self.name_matrix.append(self.name_h6)
        self.name_matrix.append(self.name_h7)
        self.name_matrix.append(self.name_h8)
        self.name_matrix.append(self.name_h9)
        self.name_matrix.append(self.name_h10)
        self.name_matrix.append(self.name_h11)
        self.name_matrix.append(self.name_h12)
        self.name_matrix.append(self.name_h13)
        self.name_matrix.append(self.name_h14)
        self.name_matrix.append(self.name_h15)
    
    
    
    
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
            '''
            self.hist_array0[i] = ROOT.TH1F(self.name_h0[i],self.name_h0[i],80,-200,200)
            self.hist_array1[i] = ROOT.TH1F(self.name_h1[i],self.name_h1[i],80,-200,200)
            self.hist_array2[i] = ROOT.TH1F(self.name_h2[i],self.name_h2[i],80,-200,200)
            self.hist_array3[i] = ROOT.TH1F(self.name_h3[i],self.name_h3[i],80,-200,200)
            self.hist_array4[i] = ROOT.TH1F(self.name_h4[i],self.name_h4[i],80,-200,200)
            self.hist_array5[i] = ROOT.TH1F(self.name_h5[i],self.name_h5[i],80,-200,200)
            self.hist_array6[i] = ROOT.TH1F(self.name_h6[i],self.name_h6[i],80,-200,200)
            self.hist_array7[i] = ROOT.TH1F(self.name_h7[i],self.name_h7[i],80,-200,200)
            self.hist_array8[i] = ROOT.TH1F(self.name_h8[i],self.name_h8[i],80,-200,200)
            self.hist_array9[i] = ROOT.TH1F(self.name_h9[i],self.name_h9[i],80,-200,200)
            self.hist_array10[i] = ROOT.TH1F(self.name_h10[i],self.name_h10[i],80,-200,200)
            self.hist_array11[i] = ROOT.TH1F(self.name_h11[i],self.name_h11[i],80,-200,200)
            self.hist_array12[i] = ROOT.TH1F(self.name_h12[i],self.name_h12[i],80,-200,200)
            self.hist_array13[i] = ROOT.TH1F(self.name_h13[i],self.name_h13[i],80,-200,200)
            self.hist_array14[i] = ROOT.TH1F(self.name_h14[i],self.name_h14[i],80,-200,200)
            self.hist_array15[i] = ROOT.TH1F(self.name_h15[i],self.name_h15[i],80,-200,200)
            '''

            self.hist_array0[i] = ROOT.TH1F(self.name_h0[i],self.name_h0[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array1[i] = ROOT.TH1F(self.name_h1[i],self.name_h1[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array2[i] = ROOT.TH1F(self.name_h2[i],self.name_h2[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array3[i] = ROOT.TH1F(self.name_h3[i],self.name_h3[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array4[i] = ROOT.TH1F(self.name_h4[i],self.name_h4[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array5[i] = ROOT.TH1F(self.name_h5[i],self.name_h5[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array6[i] = ROOT.TH1F(self.name_h6[i],self.name_h6[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array7[i] = ROOT.TH1F(self.name_h7[i],self.name_h7[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array8[i] = ROOT.TH1F(self.name_h8[i],self.name_h8[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array9[i] = ROOT.TH1F(self.name_h9[i],self.name_h9[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array10[i] = ROOT.TH1F(self.name_h10[i],self.name_h10[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array11[i] = ROOT.TH1F(self.name_h11[i],self.name_h11[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array12[i] = ROOT.TH1F(self.name_h12[i],self.name_h12[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array13[i] = ROOT.TH1F(self.name_h13[i],self.name_h13[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array14[i] = ROOT.TH1F(self.name_h14[i],self.name_h14[i],70,-TMath.Pi(),TMath.Pi())
            self.hist_array15[i] = ROOT.TH1F(self.name_h15[i],self.name_h15[i],70,-TMath.Pi(),TMath.Pi()) #80,-4,4




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
    
    def getNameMatrix(self):
        return self.name_matrix

    def getHistMatrix(self):
        return self.hist_matrix

    def getThresholds1(self):
        return self.thresholds1

    def getThresholds2(self):
        return self.thresholds2

    def getLenThresholds1(self):
        return self.lenThresholds1

    def getLenThresholds2(self):
        return self.lenThresholds2


    #def threshold_scan(self, zeroSupThreshold=5, num_events=7000):
    def threshold_scan(self, num_events=7000):
        binary_file = ixpeInputBinaryFile(FILE_PATH)
        binary_file.buildEventTable()
        numEventsInFile =  binary_file.numEvents()
        print numEventsInFile
    
        clustering = ixpeClustering(self.zeroSupThreshold)
        #clustering = ixpeClustering(5)
        #thresholds1 = range(5, 21, 1)                       # array of first pass thresholds (si ferma a 20!!!)
        #thresholds2 = range(5, 21, 1)                       # array of second pass thresholds
    
        #thresholds1 = array.array('i',(i for i in range(5,21)))
        #thresholds2 = array.array('i',(i for i in range(5,21)))
        thresholds1 = self.thresholds1
        thresholds2 = self.thresholds2
        print thresholds1[15]




        #valore = raw_input('continue?')  
    
        '''
        for i in range(0, len(thresholds1)):
            #hist_matrix.append([])
            hist_matrix[i].append(hist_array)
            #for j in range(0, 1): #len(thresholds2))
                #hist_matrix[i].append(h11)
        #valore = raw_input('continue?') 
        '''
    

        for k in range (1, numEventsInFile+1):                     # loop over all the events GIUSTO!!!!!
        #for i in range (1, num_events+1):                     # loop over all the events
            if k%1000 == 0:
                print k     
            if (k > numEventsInFile):
                break       
            evt = binary_file.readEvent(k)
            #evt = binary_file.next()
            tracks = clustering.dbScan(evt)             # ??? prende solo il cluster piu' grande? NO    # applies dbScan
            if len(tracks) != 0:            #needed! thera are events with no cluster. Is this possible? What about the trigger to have a detector event????
                track = tracks[0]                           # tracks[0] e' il cluster principale dell' i-esimo evento
                                                            # list of all pixels in the custer       # pi = pulse invariant, pha = ???

                                                  
                for i in range(0, len(thresholds1)):
                    for j in range(0,len(thresholds2)):
                        track.reconstruct(thresholds1[i], thresholds2[j], False)      # la soglia deve essere un intero (ADC)
                        phi = track.secondPassMomentsAnalysis().phi()
                        #phiDeg = ixpeMath.radToDeg(phi)
                        #self.hist_matrix[i][j].Fill(phiDeg)
                        self.hist_matrix[i][j].Fill(phi)
                        self.hist_matrix[i][j].SetTitle(self.name_matrix[i][j])
        print self.hist_matrix[15][15].GetMean()
    
        f.Write()
    
    
    
def test():
    thrScan = thresholdScan()
    thrScan.threshold_scan()
    
    
if __name__ == "__main__":

    f = TFile("threshold_scan.root", "recreate")

    #c = TCanvas("c","c",0)
    #c.cd()

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
    #print FILE_PATH

    #self.threshold_scan()

    test()    

    f.Close()


