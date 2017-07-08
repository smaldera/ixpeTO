#questo file deve essere messo nella cartella contenente le varie cartelle a loro volta contenenti le simulazioni (file recon_FITS) alle diverse energie
#usage: python -W ignore quality_cut.py sim_6keV -z 3

#M2TL e' una sorta di eccentricita' al contrario (1 = circolare, 0 = ellitticita' massima (segmento))


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
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infolder', type=str,
                        help='the input folder')
    parser.add_argument('-z', '--zero-suppression', type=int, default=5,
                        help = 'zero-suppression threshold')

    args = parser.parse_args()

    folder = args.infolder
    zero_sup = args.zero_suppression


    #ROOT file to save histograms
    f = TFile("thr_scan_%s_M2TL.root" %folder, "recreate") #NEW


    #Names and histograms matrix initialisation
    thresholds1 = array.array('i',(i for i in range(zero_sup,21)))    # thresholds are integers (ADC)
    thresholds2 = array.array('i',(i for i in range(zero_sup,21)))

    num_thresholds = 21-zero_sup

    a = np.arange(zero_sup, 21)                        # ---> per thresholds1
    b = np.arange(zero_sup, 21)                        # ---> per thresholds2
    c = np.array(['M2TL_%i,%i'%(i,j) for i,j in list(product(thresholds1,thresholds2))])
    name_matrix = np.reshape(c, (num_thresholds,num_thresholds))


    hist_array0=num_thresholds*[ROOT.TH1F]
    hist_array1=num_thresholds*[ROOT.TH1F]
    hist_array2=num_thresholds*[ROOT.TH1F]
    hist_array3=num_thresholds*[ROOT.TH1F]
    hist_array4=num_thresholds*[ROOT.TH1F]
    hist_array5=num_thresholds*[ROOT.TH1F]
    hist_array6=num_thresholds*[ROOT.TH1F]
    hist_array7=num_thresholds*[ROOT.TH1F]
    hist_array8=num_thresholds*[ROOT.TH1F]
    hist_array9=num_thresholds*[ROOT.TH1F]
    hist_array10=num_thresholds*[ROOT.TH1F]
    hist_array11=num_thresholds*[ROOT.TH1F]
    hist_array12=num_thresholds*[ROOT.TH1F]
    hist_array13=num_thresholds*[ROOT.TH1F]
    hist_array14=num_thresholds*[ROOT.TH1F]
    hist_array15=num_thresholds*[ROOT.TH1F]
    hist_array16=num_thresholds*[ROOT.TH1F]
    hist_array17=num_thresholds*[ROOT.TH1F]


    for i in range(0,num_thresholds):
        hist_array0[i] = ROOT.TH1F(name_matrix[0][i],name_matrix[0][i],70,0,1)
        hist_array1[i] = ROOT.TH1F(name_matrix[1][i],name_matrix[1][i],70,0,1)
        hist_array2[i] = ROOT.TH1F(name_matrix[2][i],name_matrix[2][i],70,0,1)
        hist_array3[i] = ROOT.TH1F(name_matrix[3][i],name_matrix[3][i],70,0,1)
        hist_array4[i] = ROOT.TH1F(name_matrix[4][i],name_matrix[4][i],70,0,1)
        hist_array5[i] = ROOT.TH1F(name_matrix[5][i],name_matrix[5][i],70,0,1)
        hist_array6[i] = ROOT.TH1F(name_matrix[6][i],name_matrix[6][i],70,0,1)
        hist_array7[i] = ROOT.TH1F(name_matrix[7][i],name_matrix[7][i],70,0,1)
        hist_array8[i] = ROOT.TH1F(name_matrix[8][i],name_matrix[8][i],70,0,1)
        hist_array9[i] = ROOT.TH1F(name_matrix[9][i],name_matrix[9][i],70,0,1)
        hist_array10[i] = ROOT.TH1F(name_matrix[10][i],name_matrix[10][i],70,0,1)
        hist_array11[i] = ROOT.TH1F(name_matrix[11][i],name_matrix[11][i],70,0,1)
        hist_array12[i] = ROOT.TH1F(name_matrix[12][i],name_matrix[12][i],70,0,1)
        hist_array13[i] = ROOT.TH1F(name_matrix[13][i],name_matrix[13][i],70,0,1)
        hist_array14[i] = ROOT.TH1F(name_matrix[14][i],name_matrix[14][i],70,0,1)
        hist_array15[i] = ROOT.TH1F(name_matrix[15][i],name_matrix[15][i],70,0,1) #80,-4,4
        hist_array16[i] = ROOT.TH1F(name_matrix[16][i],name_matrix[16][i],70,0,1)
        hist_array17[i] = ROOT.TH1F(name_matrix[17][i],name_matrix[17][i],70,0,1)


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



    for i in range (zero_sup, 21):      #thr1
        for j in range (zero_sup, 21):  #thr2
            FILE_PATH = './'+folder+'/'+folder+'_zero_%s' %zero_sup + '_th1_%s' %i + '_th2_%s' %j + '.fits'
            #print FILE_PATH

            #dataFile = open('/home/lorenzo/ixpesw/gpdsw/bin/sim_8keV/sim_8keV_L_zero_3_th1_3_th2_3.fits','r')
            dataFile = open(FILE_PATH,'r')
            fileName = dataFile.name

            #open FITS file (memmap = True to prevent RAM storage issues)
            pha_list = fits.open(fileName, memmap=True)


            #pha_list.info()
            #print(pha_list[1].columns)
        
            #load data into a separate variable
            pha_data = pha_list[1].data

    #TABLE READING TOOLS
    #  print pha_data[0]		        # row 0
    #  print pha_data[0][2]		        # row 0, element 2
    #  print pha_data[0]['track0_phi']  # row 0, element from col 'track0_phi' --> SLOWER PERFORMANCE because an intermediate row object is created
    #  print pha_data['TRK_M2L']	    # column 'track0_phi'
    #  print pha_data['track0_phi'][0]  # column 'track0_phi', element 0 --> FASTER PERFORMANCE
    #  print len(pha_data)		        # number of events

            print colored(name_matrix[i-zero_sup][j-zero_sup], 'cyan')

            num_events = len(pha_data)
            if num_events > 125000:
                num_events = 125000
            num_bins = num_events/100 #se ne metto di piu' diventa troppo lento il for loop che applica il quality cut
            #print 'num_bins = num_events/3 = ', num_bins
            #sum_num_events = 0

            hist_matrix[i-zero_sup][j-zero_sup].SetBins(num_bins,0,1)


            M2TL_array = []
            for k in range(0, num_events):
                M2TL_array.append(pha_data['TRK_M2T'][k]/pha_data['TRK_M2L'][k])

            #min_M2TL = min(M2TL_array)
            #max_M2TL = max(M2TL_array)



            #hist_M2LT = ROOT.TH1F('M2L/M2T_hist','M2L/M2T_hist',num_bins,min_M2LT,max_M2LT)
            #hist_M2LT.SetTitle('M2L/M2T_hist')
            #hist_M2LT.GetXaxis().SetTitle("M2L/M2T")
            #hist_M2LT.GetYaxis().SetTitle("counts")


            #for k in range(0, num_events):
                hist_matrix[i-zero_sup][j-zero_sup].Fill(M2TL_array[k])
                #hist_matrix[i-zero_sup][j-zero_sup].SetTitle(name_matrix[i-zero_sup][j-zero_sup])

            hist_matrix[i-zero_sup][j-zero_sup].GetXaxis().SetTitle("M2T/M2L")
            hist_matrix[i-zero_sup][j-zero_sup].GetYaxis().SetTitle("counts")
            hist_matrix[i-zero_sup][j-zero_sup].Draw()


            #valore = raw_input('continue?')


            '''
            k = 0           # A sx del primo bin (underflow bin). Il primo bin pieno e' k = 1.
            while sum_num_events < num_events*20/100:
                k = k+1
                sum_num_events = sum_num_events + hist_M2LT.GetBinContent(k)

            quality_cut = hist_M2LT.GetBinLowEdge(k+1)  # Upper edge of the last bin filled
                                                        # Lower edge of the first empty bin

            print colored('QUALITY CUT', 'green'), '= upper edge of the last filled bin (which goes over 20% of data) = '
            print '     = M2L/M2T_cut = ', quality_cut


            hist_M2LT_80 = ROOT.TH1F('M2L/M2T_hist_80','M2L/M2T_hist_80',num_bins,min_M2LT,max_M2LT)
            hist_M2LT_80.SetTitle('M2L/M2T_hist_80')
            hist_M2LT_80.GetXaxis().SetTitle("M2L/M2T")
            hist_M2LT_80.GetYaxis().SetTitle("counts")
            

            for m in range(0, num_events):
                if M2LT_array[m] >= quality_cut:
                    #fill phi histogram
                    hist_matrix[i-zero_sup][j-zero_sup].Fill(pha_data['TRK_PHI2'][m]/180*TMath.Pi())
                    hist_matrix[i-zero_sup][j-zero_sup].SetTitle(name_matrix[i-zero_sup][j-zero_sup])
                    hist_matrix[i-zero_sup][j-zero_sup].GetXaxis().SetTitle("phi (rad)")
                    hist_matrix[i-zero_sup][j-zero_sup].GetYaxis().SetTitle("counts")
                    
                    hist_M2LT_80.Fill(M2LT_array[m])
                        

            hist_matrix[i-zero_sup][j-zero_sup].Draw()
            #hist_M2LT_80.Draw()

            #print colored(name_matrix[i-zero_sup][j-zero_sup], 'cyan')
            #print colored('QUALITY CUT', 'green'), '= upper edge of the last filled bin (which goes over 20% of data) = ', quality_cut

            #valore = raw_input('continue?') 
            '''

    f.Write()    

    f.Close()

    print("--- %s seconds ---" % (time.time() - start_time))
