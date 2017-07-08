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

#usage: python -i modulation_factor_analysis.py infile.root_path -z zero_sup_thr(=3)
#python -i modulation_factor_analysis.py threshold_scan_phi.root -z zero_sup_thr(=3)


import os
import sys

import numpy
import numpy as np
import array

import matplotlib
import matplotlib.pyplot as plt

import ROOT
from ROOT import *

import math
from math import sqrt

from itertools import product

#from astropy.io import fits
#import astropy.io.fits as pyfits

from termcolor import colored

import time
start_time = time.time()


if __name__ == '__main__':


    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input phi histograms root file')
    parser.add_argument('-z', '--zero-suppression', type=int, default=3,
                        help = 'zero-suppression threshold')

    args = parser.parse_args()
    in_file = args.infile
    zero_sup    = args.zero_suppression


    # ROOT input file (phi histograms)
    inFile = TFile(in_file)


    # ROOT file to save PHI2 histograms
    outFile = TFile('mod_fact_%s' %os.path.basename(in_file).replace('thr_scan_',''), "recreate") #NEW # %.3f = float con 3 cifre dopo la virgola
    outFile.cd()


    # MATRICES INITIALIZATION
    thresholds1     = array.array('i',(i for i in range(zero_sup,21)))
    thresholds2     = array.array('i',(i for i in range(zero_sup,21)))
    len_thresholds1 = len(thresholds1)
    len_thresholds2 = len(thresholds2)
    #num_thresholds  = 21-zero_sup  #NON so se serve
    num_thresholds  = len(thresholds1)


    # Names and histograms matrix initialisation  --> PHI2
    a = np.arange(zero_sup, 21)                 # --> per thresholds1
    b = np.arange(zero_sup, 21)                 # --> per thresholds2
    c = np.array(['%i,%i'%(i,j) for i,j in list(product(thresholds1,thresholds2))])
    name_matrix = np.reshape(c, (num_thresholds,num_thresholds))



    h_modulation_factors = TH2F("h_modulation_factors", "h_modulation_factors",
    len_thresholds1, thresholds1[0]-0.5, thresholds1[len_thresholds1-1]+0.5,
    len_thresholds2, thresholds2[0]-0.5, thresholds2[len_thresholds2-1]+0.5)

    h_reduced_chi_square = TH1F("h_reduced_chi_square", "h_reduced_chi_square", 50, 0, 5)
    h_probability = TH1F("h_probability", "h_probability", 40, 0, 1) #100 bins
    h_dist_mu = TH1F("h_dist_mu", "h_dist_mu", 600, 0, 1)


    c = TCanvas("c","c",0)
    c.cd()

    maxReducedChiSquare  = 0
    maxModulationFactor  = 0
    sMaxModulationFactor = 0
    xMaxModulationFactor = 0
    yMaxModulationFactor = 0

    modulationFactors_matrix = [] #new

    for i in range (0,len_thresholds1):
        row = [] #new
        for j in range (0,len_thresholds2):
            h = inFile.Get(name_matrix[i][j])
            h.GetXaxis().SetTitle("phi (rad)")
            h.GetYaxis().SetTitle("counts")
            fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -TMath.Pi(), TMath.Pi()) #LEGGE di MALUS
            fitFunc.SetParLimits(0,0,100000000)     # offset >0
            fitFunc.SetParLimits(1,0,100000000)     # se c'e' bisogno di un'inversione, che la becchi con la fase [2]!!! 
            h.Fit("fitFunc","MR")
            gStyle.SetOptStat("nem")
            gStyle.SetOptFit(100)
            gStyle.SetStatW(0.1)
            gStyle.SetStatH(0.09)
            h.Draw("E1") #"E1" to show error bars
            h.Write() #c.Write()
            #c.Update()

            reducedChiSquare = fitFunc.GetChisquare()/fitFunc.GetNDF() #getProb
            if reducedChiSquare>maxReducedChiSquare:
                maxReducedChiSquare = reducedChiSquare

            probability = fitFunc.GetProb()

            #get covariance matrix
            fitter = TVirtualFitter.GetFitter()
            covMatrix = fitter.GetCovarianceMatrix()
            print 'cov(A,B) = cov(0,1) = ', fitter.GetCovarianceMatrixElement(0,1)
            covAB = fitter.GetCovarianceMatrixElement(0,1)


            #c.Close() #da togliere          #salvare gli istogrammi con fit su file root #quando fitto, bloccare per vedere un fit alla volta!!!
            A = fitFunc.GetParameter(0)
            B = fitFunc.GetParameter(1)
            sA = fitFunc.GetParError(0)
            sB = fitFunc.GetParError(1)
            modulationFactor = fitFunc.GetParameter(1)/(fitFunc.GetParameter(1)+2*fitFunc.GetParameter(0))
            row.append(modulationFactor) #new
            #sModulationFactor = math.sqrt( ((4*B*B*sA*sA)+(4*A*A*sB*sB)) / ((2*A+B)*(2*A+B)*(2*A+B)*(2*A+B)) )
            sModulationFactor = math.sqrt( ((4*B*B*sA*sA)+(4*A*A*sB*sB)-(8*A*B*covAB) ) / ((2*A+B)*(2*A+B)*(2*A+B)*(2*A+B)) ) #error with covariance
            print modulationFactor
            if modulationFactor>maxModulationFactor:
                maxModulationFactor  = modulationFactor
                sMaxModulationFactor = sModulationFactor
                xMaxModulationFactor = thresholds1[i]
                yMaxModulationFactor = thresholds2[j]
            h_modulation_factors.Fill(thresholds1[i],thresholds2[j],modulationFactor)
            h_reduced_chi_square.Fill(reducedChiSquare)
            h_probability.Fill(probability)

            #diffMu = 1. - modulationFactor      #Vale solo per polarizzati al 100%
                                                #Per non polarizzati, sostituire 1 con 0
            #h_dist_mu.Fill(diffMu)
            h_dist_mu.Fill(modulationFactor)


        modulationFactors_matrix.append(row) #new

    print modulationFactors_matrix[0][0] #new
    print modulationFactors_matrix[15][15] #new

    print '============================================='
    print 'Max Modulation Factor:'
    print ' 1st pass mom analysis threshold: ', xMaxModulationFactor
    print ' 2nd pass mom analysis threshold: ', yMaxModulationFactor
    print ' Max modulation factor: ', maxModulationFactor , '+-', sMaxModulationFactor
    print '============================================='

    '''
    for i in range (0,len_thresholds1):
        for j in range (0,len_thresholds2):
            if modulationFactors_matrix[i][j] > maxModulationFactor - sMaxModulationFactor:
                print '============================================='
                print 'Modulation Factors between max and max-s_max:'
                print ' 1st pass mom analysis threshold: ', thresholds1[i]
                print ' 2nd pass mom analysis threshold: ', thresholds2[j]
                print ' Modulation factor: ', modulationFactors_matrix[i][j]
                print '============================================='
                #forse devo tenere conto anche degli errori su questi???
    '''


    markerMaxModulationFactor = TMarker(xMaxModulationFactor, yMaxModulationFactor, 34)
    markerMaxModulationFactor.SetMarkerColor(2)
    markerMaxModulationFactor.SetMarkerSize(1.5)


    cc = TCanvas("cc", "cc", 0)
    cc.cd()
    h_modulation_factors.SetTitle("Modulation Factor - threshold scan")
    h_modulation_factors.GetXaxis().SetTitle("1st pass threshold")
    h_modulation_factors.GetYaxis().SetTitle("2nd pass threshold")
    h_modulation_factors.GetZaxis().SetTitle("modulation factor")
    h_modulation_factors.GetZaxis().SetTitleOffset(1.5)
    h_modulation_factors.GetZaxis().SetLabelFont(10)
    gStyle.SetOptStat(0)
    #gROOT.ForceStyle()
    #self.h_modulation_factors.UseCurrentStyle()
    #gPad.SetRightMargin(2)
    h_modulation_factors.Draw("colZ")
    markerMaxModulationFactor.Draw("same")
    cc.Write()


    ccc = TCanvas("ccc", "ccc", 0) #NEW
    ccc.cd() #NEW
    h_modulation_factors.SetTitle("Modulation Factor - threshold scan")
    gStyle.SetOptStat(0)
    #self.h_modulation_factors.UseCurrentStyle()
    h_modulation_factors.Draw("surf1Z") #NEW
    markerMaxModulationFactor.Draw("same")
    h_modulation_factors.Write()


    cccc = TCanvas("cccc", "cccc", 0)
    cccc.cd()
    h_reduced_chi_square.SetAxisRange(0, maxReducedChiSquare, "X")
    h_reduced_chi_square.GetXaxis().SetTitle("reduced chi square")
    h_reduced_chi_square.GetYaxis().SetTitle("counts")
    h_reduced_chi_square.SetTitle("Reduced Chi Square")
    #self.h_reduced_chi_square.SetBins()
    gStyle.SetStatW(0.15)
    gStyle.SetStatH(0.13)
    gStyle.SetHistFillColor(kSpring+6)
    gStyle.SetHistLineColor(kSpring-6)
    gStyle.SetHistLineStyle(3)
    gStyle.SetHistLineWidth(2)
    gStyle.SetOptStat("emr")
    h_reduced_chi_square.UseCurrentStyle()
    h_reduced_chi_square.Draw()
    h_reduced_chi_square.Write()



    ccccc = TCanvas("ccccc", "ccccc", 0)
    ccccc.cd()
    h_probability.GetXaxis().SetTitle("probability")
    h_probability.GetYaxis().SetTitle("counts")
    h_probability.SetTitle("Probability")  
    gStyle.SetStatW(0.15)
    gStyle.SetStatH(0.13)
    gStyle.SetHistFillColor(kOrange-3)
    gStyle.SetHistLineColor(kOrange+7)
    gStyle.SetHistLineStyle(3)
    gStyle.SetHistLineWidth(2)
    gStyle.SetOptStat("emr")
    h_probability.UseCurrentStyle()
    h_probability.Draw()
    h_probability.Write()


    c_dist_mu = TCanvas("cDistMu", "cDistMu", 0)
    c_dist_mu.cd()
    #h_dist_mu.GetXaxis().SetTitle("$1 - \mu$")
    h_dist_mu.GetXaxis().SetTitle("$\mu$")
    h_dist_mu.GetYaxis().SetTitle("counts")
    #h_dist_mu.SetTitle("$\mu_{real} - \mu$")  
    h_dist_mu.SetTitle("$\mu$ distribution")
    gStyle.SetStatW(0.15)
    gStyle.SetStatH(0.13)
    gStyle.SetHistFillColor(kBlue-4)
    gStyle.SetHistLineColor(kBlue+3)
    gStyle.SetHistLineStyle(3)
    gStyle.SetHistLineWidth(2)
    gStyle.SetOptStat("emr")
    gStyle.SetStatW(0.16)
    gStyle.SetStatH(0.12)
    h_dist_mu.UseCurrentStyle()
    h_dist_mu.Draw()
    h_dist_mu.Write()

    #valore = raw_input('continue?')
    
    outFile.Close()
    inFile.Close()
    



    print("--- %s seconds ---" % (time.time() - start_time))

        #TO DO:

        #draw TH2 in 3-d??? (surf2) DONE
        #set Stat per chi quadro (opzione di TH1/2 o canvas?? altrimenti da TBrowser View-->Editor-->Chi!! DONE
        #istogramma del chi quadro ridotto!! DONE
        #che errore da' di default sull'altezza delle barre dell'istogramma (sembra sqrt(N_in_bin))
        #soglie t.c. mod_factor == max ---> draw on the TH2 --> write canvas con h_modulation_factors + TMarker DONE
        #angolo fisso, diverse energie DONE
        #energia fissa, diversi angoli
        #con diverso asse di polarizzazione non dovrebbe cambiare
        #fissa range di variazione colore sul TH2 NO!
        #per vedere la dipendenza dalle due soglie (quanto dipende dalla variazione della prima? quanto della seconda?) DONE
        #       --> possiamo sommare tutti i valori del modulation_factor corrismondenti alla stessa thr1 e per le diverse thr2 DONE
        #sistemare nomi file output (DONE) e i titoli dei grafici 


        #ERRORE dul fattore di modulazione?
        #e se il massimo assoluto segnato con il TMarker fosse una fluttuazione? Lontano dalla zona di crescita del mu?
        #andamento del chi quadro in funzione delle soglie
        #il massimo del fattore di modulazione e' nella zona dove il chi quadro ridotto e' buono?
        # vediamo...

        #potrei vedere anche come varia il punto di conversione. rimane stabile li' vicino o si sposta???
        #queste due soglie cosa fanno effettivamente nell'algoritmo?

        #cambiare i nomi degli istogrammi (name_matrix) --> non soglia numero 1,1  , bensi' 5,5 (valore in ADC delle soglie!) DONE

        #METTERE DEI TAGLI SUL CHI QUADRO O SULLA PROBABILITA'???    



