#questo file deve essere messo nella cartella contenente thr_scan_sim_6keV_M2TL.root (alle diverse energie) e le cartelle con i vari FITS file alle diverse energie
#usage: python -W ignore set_quality_cut.py sim_6keV path/to/M2TL_root_file.root -z 3


import os
import sys

import numpy as np
import array
from itertools import product

import matplotlib
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.io.fits as pyfits

import ROOT
from ROOT import *
from math import sqrt

from termcolor import colored

import time
start_time = time.time()


if __name__ == '__main__':

    import argparse
    formatter   = argparse.ArgumentDefaultsHelpFormatter
    parser      = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infolder', type=str, help='the input folder')
    parser.add_argument('infile', type=str, help = 'the input M2TL root file')
    parser.add_argument('-z', '--zero-suppression', type=int, default=5, help = 'zero-suppression threshold')

    args        = parser.parse_args()
    folder      = args.infolder
    in_file     = args.infile
    zero_sup    = args.zero_suppression


    # ROOT input file (M2TL histograms)
    inFile_M2TL = TFile(in_file)


    # ROOT file to save PHI2 histograms
    f = TFile("thr_scan_%s_PHI2.root" %folder, "recreate") #NEW


    # MATRICES INITIALIZATION
    thresholds1     = array.array('i',(i for i in range(zero_sup,21)))
    thresholds2     = array.array('i',(i for i in range(zero_sup,21)))
    num_thresholds  = 21-zero_sup


    # Names and histograms matrix initialisation  --> M2TL (from root file)
    aM2TL   = np.arange(zero_sup, 21)           # --> per thresholds1
    bM2TL   = np.arange(zero_sup, 21)           # --> per thresholds2
    cM2TL   = np.array(['M2TL_%i,%i'%(i,j) for i,j in list(product(thresholds1,thresholds2))])
    name_matrix_M2TL = np.reshape(cM2TL, (num_thresholds,num_thresholds))


    # Names and histograms matrix initialisation  --> PHI2
    a = np.arange(zero_sup, 21)                 # --> per thresholds1
    b = np.arange(zero_sup, 21)                 # --> per thresholds2
    c = np.array(['%i,%i'%(i,j) for i,j in list(product(thresholds1,thresholds2))])
    name_matrix = np.reshape(c, (num_thresholds,num_thresholds))


    hist_array0 =num_thresholds*[ROOT.TH1F]
    hist_array1 =num_thresholds*[ROOT.TH1F]
    hist_array2 =num_thresholds*[ROOT.TH1F]
    hist_array3 =num_thresholds*[ROOT.TH1F]
    hist_array4 =num_thresholds*[ROOT.TH1F]
    hist_array5 =num_thresholds*[ROOT.TH1F]
    hist_array6 =num_thresholds*[ROOT.TH1F]
    hist_array7 =num_thresholds*[ROOT.TH1F]
    hist_array8 =num_thresholds*[ROOT.TH1F]
    hist_array9 =num_thresholds*[ROOT.TH1F]
    hist_array10=num_thresholds*[ROOT.TH1F]
    hist_array11=num_thresholds*[ROOT.TH1F]
    hist_array12=num_thresholds*[ROOT.TH1F]
    hist_array13=num_thresholds*[ROOT.TH1F]
    hist_array14=num_thresholds*[ROOT.TH1F]
    hist_array15=num_thresholds*[ROOT.TH1F]
    hist_array16=num_thresholds*[ROOT.TH1F]
    hist_array17=num_thresholds*[ROOT.TH1F]


    for i in range(0, num_thresholds):
        hist_array0[i]  = ROOT.TH1F(name_matrix[0][i], name_matrix[0][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array1[i]  = ROOT.TH1F(name_matrix[1][i], name_matrix[1][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array2[i]  = ROOT.TH1F(name_matrix[2][i], name_matrix[2][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array3[i]  = ROOT.TH1F(name_matrix[3][i], name_matrix[3][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array4[i]  = ROOT.TH1F(name_matrix[4][i], name_matrix[4][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array5[i]  = ROOT.TH1F(name_matrix[5][i], name_matrix[5][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array6[i]  = ROOT.TH1F(name_matrix[6][i], name_matrix[6][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array7[i]  = ROOT.TH1F(name_matrix[7][i], name_matrix[7][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array8[i]  = ROOT.TH1F(name_matrix[8][i], name_matrix[8][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array9[i]  = ROOT.TH1F(name_matrix[9][i], name_matrix[9][i], 70,-TMath.Pi(),TMath.Pi())
        hist_array10[i] = ROOT.TH1F(name_matrix[10][i],name_matrix[10][i],70,-TMath.Pi(),TMath.Pi())
        hist_array11[i] = ROOT.TH1F(name_matrix[11][i],name_matrix[11][i],70,-TMath.Pi(),TMath.Pi())
        hist_array12[i] = ROOT.TH1F(name_matrix[12][i],name_matrix[12][i],70,-TMath.Pi(),TMath.Pi())
        hist_array13[i] = ROOT.TH1F(name_matrix[13][i],name_matrix[13][i],70,-TMath.Pi(),TMath.Pi())
        hist_array14[i] = ROOT.TH1F(name_matrix[14][i],name_matrix[14][i],70,-TMath.Pi(),TMath.Pi())
        hist_array15[i] = ROOT.TH1F(name_matrix[15][i],name_matrix[15][i],70,-TMath.Pi(),TMath.Pi())
        hist_array16[i] = ROOT.TH1F(name_matrix[16][i],name_matrix[16][i],70,-TMath.Pi(),TMath.Pi())
        hist_array17[i] = ROOT.TH1F(name_matrix[17][i],name_matrix[17][i],70,-TMath.Pi(),TMath.Pi())


    hist_matrix = []

    hist_matrix.append(hist_array0)
    hist_matrix.append(hist_array1)
    hist_matrix.append(hist_array2)
    hist_matrix.append(hist_array3)
    hist_matrix.append(hist_array4)
    hist_matrix.append(hist_array5)
    hist_matrix.append(hist_array6)
    hist_matrix.append(hist_array7)
    hist_matrix.append(hist_array8)
    hist_matrix.append(hist_array9)
    hist_matrix.append(hist_array10)
    hist_matrix.append(hist_array11)
    hist_matrix.append(hist_array12)
    hist_matrix.append(hist_array13)
    hist_matrix.append(hist_array14)
    hist_matrix.append(hist_array15)
    hist_matrix.append(hist_array16)
    hist_matrix.append(hist_array17)



    # Loop over thresholds couples
    for i in range (zero_sup, 21):          #thr1
        for j in range (zero_sup, 21):      #thr2
            print colored(name_matrix[i-zero_sup][j-zero_sup], 'cyan')

            FITS_FILE_PATH = './'+folder+'/'+folder+'_zero_%s' %zero_sup + '_th1_%s' %i + '_th2_%s' %j + '.fits'

            # Open FITS file
            dataFile = open(FITS_FILE_PATH,'r')
            fileName = dataFile.name
            pha_list = fits.open(fileName, memmap=True)
            pha_data = pha_list[1].data

            num_events = len(pha_data)
            if num_events > 125000:
                num_events = 125000
            num_bins = num_events/100
            sum_num_events = 0


            # Get M2TL histogram from ROOT M2TL input file
            #inFile_M2TL = TFile(in_file) --> togliere
            hist_M2TL = inFile_M2TL.Get(name_matrix_M2TL[i-zero_sup][j-zero_sup])

            
            # Quality cut
            k = num_bins+1           # A dx dell'ultimo bin (overflow bin). L'ultimo bin pieno e' k = nbins.
            while sum_num_events < num_events*20/100:
                k = k-1
                sum_num_events = sum_num_events + hist_M2TL.GetBinContent(k)

            quality_cut = hist_M2TL.GetBinLowEdge(k)

            print colored('QUALITY CUT', 'green'), '= lower edge of the last filled bin (which goes over 20% of data) = '
            print '     = M2L/M2T_cut = ', quality_cut


            #hist_M2TL_80 = ROOT.TH1F('M2T/M2L_hist_80','M2T/M2L_hist_80',num_bins,0,1)
            #hist_M2TL_80.SetTitle('M2T/M2L_hist_80')
            #hist_M2TL_80.GetXaxis().SetTitle("M2T/M2L")
            #hist_M2TL_80.GetYaxis().SetTitle("counts")


            # Phi histograms (only events which pass the quality cut)
            M2TL_array = []
            for m in range(0, num_events):
                M2TL_array.append(pha_data['TRK_M2T'][m]/pha_data['TRK_M2L'][m])
                if M2TL_array[m] <= quality_cut:
                    # Fill PHI histogram
                    hist_matrix[i-zero_sup][j-zero_sup].Fill(pha_data['TRK_PHI2'][m]/180*TMath.Pi())
                    hist_matrix[i-zero_sup][j-zero_sup].SetTitle(name_matrix[i-zero_sup][j-zero_sup])
                    hist_matrix[i-zero_sup][j-zero_sup].GetXaxis().SetTitle("TRK_PHI2 (rad)") # TRK_PHI2 e' il phi del python wrap
                    hist_matrix[i-zero_sup][j-zero_sup].GetYaxis().SetTitle("counts")
                    
                    #hist_M2TL_80.Fill(M2TL_array[m])
                        

            hist_matrix[i-zero_sup][j-zero_sup].Draw()
            #hist_M2TL_80.Draw()
            
            #print("--- %s seconds ---" % (time.time() - start_time))
            #valore = raw_input('continue?') 


    f.Write()    

    f.Close()

    print("--- %s seconds ---" % (time.time() - start_time))




    # FITS TABLES reading tools
    #  print pha_data['TRK_M2L']	    # column 'track0_phi'
    #  print pha_data['track0_phi'][0]  # column 'track0_phi', element 0 --> FASTER PERFORMANCE
    #  print len(pha_data)		        # number of events

