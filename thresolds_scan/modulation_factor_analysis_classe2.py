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

        self.len_thresholds1 = thresholdScanObj.lenThresholds1
        self.len_thresholds2 = thresholdScanObj.lenThresholds2
        self.name_matrix = thresholdScanObj.name_matrix

        #self.c_init=0
        self.h_modulation_factors = TH2F()
        self.h_reduced_chi_square = TH1F()
        self.markerMaxModulationFactor = TMarker()

    def histogram_fitter(self):
        

        name_matrix = self.name_matrix
        thresholds1 = self.thresholds1
        thresholds2 = self.thresholds2
        len_thresholds1 = self.len_thresholds1
        len_thresholds2 = self.len_thresholds2

    
        inFile = TFile(FILE_PATH)

        #self.c_init.Clear()
        #self.c_init.cd()

        self.h_modulation_factors = TH2F("h_modulation_factors", "h_modulation_factors", len_thresholds1, thresholds1[0]-0.5, thresholds1[len_thresholds1-1]+0.5, len_thresholds2, thresholds2[0]-0.5, thresholds2[len_thresholds2-1]+0.5)
        self.h_reduced_chi_square = TH1F("h_reduced_chi_square", "h_reduced_chi_square", 50, 0, 5)

        out_file = TFile('mod_fact_%s' %os.path.basename(FILE_PATH).replace('thr_scan_',''), "recreate") #NEW # %.3f = float con 3 cifre dopo la virgola
        #out_file = TFile("out_file.root", "recreate") #OLD
        out_file.cd()

        c = TCanvas("c","c",0)
        c.cd()

        maxReducedChiSquare  = 0
        maxModulationFactor  = 0
        xMaxModulationFactor = 0
        yMaxModulationFactor = 0

        for i in range (0,len_thresholds1):
            for j in range (0,len_thresholds2):

                h = inFile.Get(name_matrix[i][j])
                fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -TMath.Pi(), TMath.Pi())
                fitFunc.SetParLimits(0,0,100000000)     # offset >0
                fitFunc.SetParLimits(1,0,100000000)     # se c'e' bisogno di un'inversione, che la becchi con la fase [2]!!! 
                h.Fit("fitFunc","MR")
                gStyle.SetOptStat("nem")
                gStyle.SetOptFit(100)
                gStyle.SetStatW(0.1)
                gStyle.SetStatH(0.09)
                h.Draw("E1") #"E1" to show error bars
                c.Update()
                h.Write()
                reducedChiSquare = fitFunc.GetChisquare()/fitFunc.GetNDF()
                if reducedChiSquare>maxReducedChiSquare:
                    maxReducedChiSquare = reducedChiSquare

                #c.Close() #da togliere          #salvare gli istogrammi con fit su file root #quando fitto, bloccare per vedere un fit alla volta!!!
                modulationFactor = fitFunc.GetParameter(1)/(fitFunc.GetParameter(1)+2*fitFunc.GetParameter(0)) #NO!!! Da cambiare!!!!
                print modulationFactor
                if modulationFactor>maxModulationFactor:
                    maxModulationFactor  = modulationFactor
                    xMaxModulationFactor = thresholds1[i]
                    yMaxModulationFactor = thresholds2[j]
                self.h_modulation_factors.Fill(thresholds1[i],thresholds2[j],modulationFactor)
                self.h_reduced_chi_square.Fill(reducedChiSquare)

        self.markerMaxModulationFactor = TMarker(xMaxModulationFactor, yMaxModulationFactor, 34)
        self.markerMaxModulationFactor.SetMarkerColor(2)
        self.markerMaxModulationFactor.SetMarkerSize(1.5)

        cc = TCanvas("cc", "cc", 0)
        cc.cd()
        self.h_modulation_factors.SetTitle("Modulation Factor - threshold scan")
        self.h_modulation_factors.GetXaxis().SetTitle("1st pass threshold")
        self.h_modulation_factors.GetYaxis().SetTitle("2nd pass threshold")
        self.h_modulation_factors.GetZaxis().SetTitle("modulation factor")
        self.h_modulation_factors.GetZaxis().SetTitleOffset(1.5)
        self.h_modulation_factors.GetZaxis().SetLabelFont(10)
        gStyle.SetOptStat(0)
        #gROOT.ForceStyle()
        #self.h_modulation_factors.UseCurrentStyle()
        #gPad.SetRightMargin(2)
        self.h_modulation_factors.Draw("colZ")
        self.markerMaxModulationFactor.Draw("same")
        cc.Write()

        ccc = TCanvas("ccc", "ccc", 0) #NEW
        ccc.cd() #NEW
        self.h_modulation_factors.SetTitle("Modulation Factor - threshold scan")
        gStyle.SetOptStat(0)
        #self.h_modulation_factors.UseCurrentStyle()
        self.h_modulation_factors.Draw("surf1Z") #NEW
        self.markerMaxModulationFactor.Draw("same")
        self.h_modulation_factors.Write()

        cccc = TCanvas("cccc", "cccc", 0)
        cccc.cd()
        self.h_reduced_chi_square.SetAxisRange(0, maxReducedChiSquare, "X")
        self.h_reduced_chi_square.SetTitle("Reduced Chi Square")
        #self.h_reduced_chi_square.SetBins()
        gStyle.SetStatW(0.15)
        gStyle.SetStatH(0.13)
        gStyle.SetHistFillColor(kSpring+6)
        gStyle.SetHistLineColor(kSpring-6)
        gStyle.SetHistLineStyle(3)
        gStyle.SetHistLineWidth(2)
        gStyle.SetOptStat("emr")
        self.h_reduced_chi_square.UseCurrentStyle()
        self.h_reduced_chi_square.Draw()
        self.h_reduced_chi_square.Write()
    
        valore = raw_input('continue?')
    
        out_file.Close()
        inFile.Close()
    
        #TO DO:

        #draw TH2 in 3-d??? (surf2) DONE
        #set Stat per chi quadro (opzione di TH1/2 o canvas?? altrimenti da TBrowser View-->Editor-->Chi!! DONE
        #istogramma del chi quadro ridotto!! 
        #che errore da' di default sull'altezza delle barre dell'istogramma (sembra sqrt(N_in_bin))
        #soglie t.c. mod_factor == max ---> draw on the TH2 --> write canvas con h_modulation_factors + TMarker DONE
        #angolo fisso, diverse energie
        #energia fissa, diversi angoli
        #con diverso asse di polarizzazione non dovrebbe cambiare
        #fissa range di variazione colore sul TH2
        #per vedere la dipendenza dalle due soglie (quanto dipende dalla variazione della prima? quanto della seconda?)
        #       --> possiamo sommare tutti i valori del modulation_factor corrismondenti alla stessa thr1 e per le diverse thr2
        #sistemare nomi file output (DONE) e i titoli dei grafici



def test():

    #c = TCanvas("c","c",0)

    thresholdScanObj = thresholdScan()  # oggetto della classe thresholdScan
                                        # definita in mom_analysis_threshold_scan_classe.py

    modFact = modulationFactor(thresholdScanObj)    # oggetto della classe modulationFactor
                                                    # definita in questo file
    #modFact.c_init = c
    modFact.histogram_fitter()



if __name__ == "__main__":

    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input root file')

    args = parser.parse_args()

    FILE_PATH = args.infile

    test()
    



