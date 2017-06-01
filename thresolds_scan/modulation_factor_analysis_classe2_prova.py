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

#usage: python -i modulation_factor_analysis_classe2.py infile.root_path
#python -i modulation_factor_analysis_classe2.py threshold_scan.root


import os
import numpy
import array

import matplotlib.pyplot as plt

import ROOT
from ROOT import *

from gpdswswig.Recon import *
from gpdswswig.Utils import ixpeMath
from gpdswswig.Io import ixpeInputBinaryFile

import mom_analysis_threshold_scan_classe as mom
from mom_analysis_threshold_scan_classe import *



class modulationFactor:

    def __init__(self, thresholdScanObj):

        self.thresholds1 = thresholdScanObj.thresholds1
        self.thresholds2 = thresholdScanObj.thresholds2
        '''
        self.len_thresholds1 = thresholdScanObj.getLenThresholds1()
        self.len_thresholds2 = thresholdScanObj.getLenThresholds2()
        '''
        self.len_thresholds1 = thresholdScanObj.lenThresholds1
        self.len_thresholds2 = thresholdScanObj.lenThresholds2
        self.name_matrix = thresholdScanObj.name_matrix
        #self.hist_matrix = thresholdScanObj.hist_matrix

        self.c_init=0
        self.h_modulation_factors = TH2F()




        '''
        self.hist0 = [h1_1,h1_2,h1_3,h1_4,h1_5,h1_6,h1_7,h1_8,h1_9,h1_10,h1_11,h1_12,h1_13,h1_14,h1_15,h1_16]
        self.hist1 = [h2_1,h2_2,h2_3,h2_4,h2_5,h2_6,h2_7,h2_8,h2_9,h2_10,h2_11,h2_12,h2_13,h2_14,h2_15,h2_16]
        self.hist2 = [h3_1,h3_2,h3_3,h3_4,h3_5,h3_6,h3_7,h3_8,h3_9,h3_10,h3_11,h3_12,h3_13,h3_14,h3_15,h3_16]
        self.hist3 = [h4_1,h4_2,h4_3,h4_4,h4_5,h4_6,h4_7,h4_8,h4_9,h4_10,h4_11,h4_12,h4_13,h4_14,h4_15,h4_16]
        self.hist4 = [h5_1,h5_2,h5_3,h5_4,h5_5,h5_6,h5_7,h5_8,h5_9,h5_10,h5_11,h5_12,h5_13,h5_14,h5_15,h5_16]
        self.hist5 = [h6_1,h6_2,h6_3,h6_4,h6_5,h6_6,h6_7,h6_8,h6_9,h6_10,h6_11,h6_12,h6_13,h6_14,h6_15,h6_16]
        self.hist6 = [h7_1,h7_2,h7_3,h7_4,h7_5,h7_6,h7_7,h7_8,h7_9,h7_10,h7_11,h7_12,h7_13,h7_14,h7_15,h7_16]
        self.hist7 = [h8_1,h8_2,h8_3,h8_4,h8_5,h8_6,h8_7,h8_8,h8_9,h8_10,h8_11,h8_12,h8_13,h8_14,h8_15,h8_16]
        self.hist8 = [h9_1,h9_2,h9_3,h9_4,h9_5,h9_6,h9_7,h9_8,h9_9,h9_10,h9_11,h9_12,h9_13,h9_14,h9_15,h9_16]
        self.hist9 = [h10_1,h10_2,h10_3,h10_4,h10_5,h10_6,h10_7,h10_8,h10_9,h10_10,h10_11,h10_12,h10_13,h10_14,h10_15,h10_16]
        self.hist10= [h11_1,h11_2,h11_3,h11_4,h11_5,h11_6,h11_7,h11_8,h11_9,h11_10,h11_11,h11_12,h11_13,h11_14,h11_15,h11_16]
        self.hist11= [h12_1,h12_2,h12_3,h12_4,h12_5,h12_6,h12_7,h12_8,h12_9,h12_10,h12_11,h12_12,h12_13,h12_14,h12_15,h12_16]
        self.hist12= [h13_1,h13_2,h13_3,h13_4,h13_5,h13_6,h13_7,h13_8,h13_9,h13_10,h13_11,h13_12,h13_13,h13_14,h13_15,h13_16]
        self.hist13= [h14_1,h14_2,h14_3,h14_4,h14_5,h14_6,h14_7,h14_8,h14_9,h14_10,h14_11,h14_12,h14_13,h14_14,h14_15,h14_16]
        self.hist14= [h15_1,h15_2,h15_3,h15_4,h15_5,h15_6,h15_7,h15_8,h15_9,h15_10,h15_11,h15_12,h15_13,h15_14,h15_15,h15_16]
        self.hist15= [h16_1,h16_2,h16_3,h16_4,h16_5,h16_6,h16_7,h16_8,h16_9,h16_10,h16_11,h16_12,h16_13,h16_14,h16_15,h16_16]
    
        self.h_matrix = []
        self.h_matrix.append(self.hist0)
        self.h_matrix.append(self.hist1)
        self.h_matrix.append(self.hist2)
        self.h_matrix.append(self.hist3)
        self.h_matrix.append(self.hist4)
        self.h_matrix.append(self.hist5)
        self.h_matrix.append(self.hist6)
        self.h_matrix.append(self.hist7)
        self.h_matrix.append(self.hist8)
        self.h_matrix.append(self.hist9)
        self.h_matrix.append(self.hist10)
        self.h_matrix.append(self.hist11)
        self.h_matrix.append(self.hist12)
        self.h_matrix.append(self.hist13)
        self.h_matrix.append(self.hist14)
        self.h_matrix.append(self.hist15)
        '''


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







    def histogram_fitter(self):
        
        #valore = raw_input('continue?')

        #momAnalysisThresholdScan = thresholdScan() #from mom_analysis_threshold_scan_class.py
        #name_matrix = momAnalysisThresholdScan.getNameMatrix()
        name_matrix = self.name_matrix
        hist_matrix = self.hist_matrix
        thresholds1 = self.thresholds1
        thresholds2 = self.thresholds2
        len_thresholds1 = self.len_thresholds1
        len_thresholds2 = self.len_thresholds2

    
        inFile = TFile(FILE_PATH)
        #print inFile.ls()
    
        #GET hist from input file
        self.c_init.Clear()
        #cc = TCanvas("cc","cc",0)
        self.c_init.cd()

        self.h_modulation_factors = TH2F("h_modulation_factors", "h_modulation_factors", len_thresholds1, thresholds1[0]-0.5, thresholds1[len_thresholds1-1]+0.5, len_thresholds2, thresholds2[0]-0.5, thresholds2[len_thresholds2-1]+0.5)

    
        for i in range (0,len_thresholds1):
            for j in range (0,len_thresholds2):
                #DA PROVARE: usare self.hist_matrix dall'altra classe, cosi' sovrascrivo gli istogrammi e posso salvarli su ROOT
                #h = inFile.Get(name_matrix[i][j])
                hist_matrix[i][j] = inFile.Get(name_matrix[i][j])
                fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -4, 4)
                #fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -180, 180)  #LBaldini, presentazione per csn2, p. 11
                                                                                          #e' giusto inserire cosi' la dilatazione [3]?????
                                                                                          #sballa il conto per il fattore di modulazione che c'e' sulla slide
                #fitFunc.SetParLimits(0,0,100000000)
                fitFunc.SetParLimits(1,0,100000000)           #0 > 0, 1 > 0
                c = TCanvas("c","c",0) #si puo' pensare di definirlo altrove
                c.cd()
                hist_matrix[i][j].Fit("fitFunc","MR")   # opzione N per non disegnare
                #valore = raw_input('continue?')
                #c.Close() #da togliere          #salvare gli istogrammi con fit su file root #quando fitto, bloccare per vedere un fit alla volta!!!
                modulationFactor = fitFunc.GetParameter(1)/(fitFunc.GetParameter(1)+2*fitFunc.GetParameter(0)) #NO!!! Da cambiare!!!!
                print modulationFactor
                cc.cd()
                h_modulation_factors.Fill(thresholds1[i],thresholds2[j],modulationFactor)
                h_modulation_factors.GetXaxis().SetTitle("1st pass threshold")
                h_modulation_factors.GetYaxis().SetTitle("2nd pass threshold")
                h_modulation_factors.Draw("colZ")
    
        valore = raw_input('continue?')
    
        ff.Write()

        inFile.Close()
    


def test():

    cc = TCanvas("cc","cc",0)

    thresholdScanObj = thresholdScan()  # oggetto della classe thresholdScan
                                        # definita in mom_analysis_threshold_scan_classe.py

    modFact = modulationFactor(thresholdScanObj)    # oggetto della classe modulationFactor
                                                    # definita in questo file
    modFact.c_init = cc
    modFact.histogram_fitter()



if __name__ == "__main__":

    ff = TFile("modulation_factor.root", "recreate")

    #c = TCanvas("c","c",0)
    #c.cd()

    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input root file')


    args = parser.parse_args()

    FILE_PATH = args.infile


    test()
    
    ff.Close()


