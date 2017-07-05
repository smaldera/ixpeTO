# script per leggere i file .fits
# usage: python -i read_fits.py path/to/FITSfile.fits output_file_name.root


#questo file deve essere messo nella cartella contenente le varie cartelle a loro volta contenenti le simulazioni (file recon_FITS) alle diverse energie
#usage: python quality_cut.py sim_6keV -z 3


import os
import sys

import numpy as np
import array
#from array import array
from itertools import product

import matplotlib
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.io.fits as pyfits

import ROOT
from ROOT import *
from math import sqrt

#from gpdswswig.Recon import *
#from gpdswswig.Utils import ixpeMath
#from gpdswswig.Io import ixpeInputBinaryFile



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
    
    #print os.path.basename(FILE_PATH)
    #print 'ZERO_SUPPRESSION set to ', zero_sup, ' ADC'
    print folder

    #cmd = 'cd ' + folder
    #os.system(cmd)




    f = TFile("thr_scan_%s_PHI2.root" %folder, "recreate") #NEW






    thresholds1 = array.array('i',(i for i in range(zero_sup,21)))    # thresholds are integers (ADC)
    thresholds2 = array.array('i',(i for i in range(zero_sup,21)))

    num_thresholds = 21-zero_sup

    #Initialize name matrix
    a = np.arange(zero_sup, 21)                        # ---> per thresholds1
    b = np.arange(zero_sup, 21)                        # ---> per thresholds2
    c = np.array(['%i,%i'%(i,j) for i,j in list(product(thresholds1,thresholds2))])
    name_matrix = np.reshape(c, (num_thresholds,num_thresholds))


    #Initialize histograms matrix
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
        hist_array0[i] = ROOT.TH1F(name_matrix[0][i],name_matrix[0][i],70,-TMath.Pi(),TMath.Pi())
        hist_array1[i] = ROOT.TH1F(name_matrix[1][i],name_matrix[1][i],70,-TMath.Pi(),TMath.Pi())
        hist_array2[i] = ROOT.TH1F(name_matrix[2][i],name_matrix[2][i],70,-TMath.Pi(),TMath.Pi())
        hist_array3[i] = ROOT.TH1F(name_matrix[3][i],name_matrix[3][i],70,-TMath.Pi(),TMath.Pi())
        hist_array4[i] = ROOT.TH1F(name_matrix[4][i],name_matrix[4][i],70,-TMath.Pi(),TMath.Pi())
        hist_array5[i] = ROOT.TH1F(name_matrix[5][i],name_matrix[5][i],70,-TMath.Pi(),TMath.Pi())
        hist_array6[i] = ROOT.TH1F(name_matrix[6][i],name_matrix[6][i],70,-TMath.Pi(),TMath.Pi())
        hist_array7[i] = ROOT.TH1F(name_matrix[7][i],name_matrix[7][i],70,-TMath.Pi(),TMath.Pi())
        hist_array8[i] = ROOT.TH1F(name_matrix[8][i],name_matrix[8][i],70,-TMath.Pi(),TMath.Pi())
        hist_array9[i] = ROOT.TH1F(name_matrix[9][i],name_matrix[9][i],70,-TMath.Pi(),TMath.Pi())
        hist_array10[i] = ROOT.TH1F(name_matrix[10][i],name_matrix[10][i],70,-TMath.Pi(),TMath.Pi())
        hist_array11[i] = ROOT.TH1F(name_matrix[11][i],name_matrix[11][i],70,-TMath.Pi(),TMath.Pi())
        hist_array12[i] = ROOT.TH1F(name_matrix[12][i],name_matrix[12][i],70,-TMath.Pi(),TMath.Pi())
        hist_array13[i] = ROOT.TH1F(name_matrix[13][i],name_matrix[13][i],70,-TMath.Pi(),TMath.Pi())
        hist_array14[i] = ROOT.TH1F(name_matrix[14][i],name_matrix[14][i],70,-TMath.Pi(),TMath.Pi())
        hist_array15[i] = ROOT.TH1F(name_matrix[15][i],name_matrix[15][i],70,-TMath.Pi(),TMath.Pi()) #80,-4,4
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


    '''    
    for i in range (zero_sup, 21):
        for j in range (zero_sup, 21):
            cmd='./ixperecon --input-files ' + FILE_PATH + ' --output-folder .' + ' --output-suffix zero_%s' %zero_sup + '_th1_%s' %i + '_th2_%s' %j + ' --threshold %s' %zero_sup + ' --moma1-threshold %s' %i + ' --moma2-threshold %s' %j #+str(i)
            print "sto per eseguire ",cmd
            os.system(cmd)


    #if number of args not correct
      if (len(sys.argv) != 3): #2
        print "    usage: python -i read_fits.py path/to/FITSfile.fits output_file_name.root"
        sys.exit()
    
    #  print sys.argv[2]
    '''



    #get FITS file name
    #FILE_PATH = '/'+folder+'/'+folder+'_L_zero_%s' %zero_sup + '_th1_%s' %3 + '_th2_%s' %3 + '.fits'
    for i in range (zero_sup, 21): #thr1
        for j in range (zero_sup, 21): #thr2
            FILE_PATH = './'+folder+'/'+folder+'_L_zero_%s' %zero_sup + '_th1_%s' %i + '_th2_%s' %j + '.fits'
            print FILE_PATH

            #dataFile = open('/home/lorenzo/ixpesw/gpdsw/bin/sim_8keV/sim_8keV_L_zero_3_th1_3_th2_3.fits','r')
            #dataFile = open('./sim_8keV/sim_8keV_L_zero_3_th1_3_th2_3.fits','r')
            dataFile = open(FILE_PATH,'r')
            #dataFile = open(sys.argv[1],'r')
            fileName = dataFile.name

            #open FITS file (memmap = True to prevent RAM storage issues)
            pha_list = fits.open(fileName, memmap=True)

            #print info on FITS file content
            pha_list.info()
    
            #print columns name, format, unit
            print(pha_list[1].columns)
        
            #load data into a separate variable
            pha_data = pha_list[1].data

            #valore = raw_input('continue?')    

    #TABLE READING TOOLS
    #  print pha_data[0]		   # row 0
    #  print pha_data[0][2]		   # row 0, element 2
    #  print pha_data[0]['track0_phi']  # row 0, element from col 'track0_phi' --> SLOWER PERFORMANCE because an intermediate row object gets created
            #print pha_data['TRK_M2L']	   # column 'track0_phi'
    #  print pha_data['track0_phi'][0]  # column 'track0_phi', element 0 --> FASTER PERFORMANCE
    #  print len(pha_data)		   # number of events


            num_events = len(pha_data)
            num_bins = num_events/3


            sum_num_events = 0


            M2LT_array = []
            for k in range(0, num_events):
                M2LT_array.append(pha_data['TRK_M2L'][k]/pha_data['TRK_M2T'][k])

            min_M2LT = min(M2LT_array)
            max_M2LT = max(M2LT_array)

            hist_M2LT = ROOT.TH1F('M2L/M2T_hist','M2L/M2T_hist',num_bins,min_M2LT,max_M2LT)
            hist_M2LT.SetTitle('M2L/M2T_hist')
            hist_M2LT.GetXaxis().SetTitle("M2L/M2T")
            hist_M2LT.GetYaxis().SetTitle("counts")

            #if sum_num_events < num_events*20/100:

            '''
            k = 0
            while sum_num_events < 12000: #num_events*20/100:
                hist_M2LT.Fill(M2LT_array[k])
                sum_num_events = sum_num_events + hist_M2LT.GetBinContent(k+1)
                k = k+1

            hist_M2LT.Draw()
            '''


            for k in range(0, num_events):
                hist_M2LT.Fill(M2LT_array[k])

            hist_M2LT.Draw()

            k = 0 #a sinistra del primo bin
            while sum_num_events < num_events*20/100:
                k = k+1
                sum_num_events = sum_num_events + hist_M2LT.GetBinContent(k)


            quality_cut = hist_M2LT.GetBinLowEdge(k+1)

            print 'Upper edge of the bin which goes over 20% = quality cut: ', quality_cut


            for m in range(0, num_events):
                if M2LT_array[m] >= quality_cut:
                    #fill phi histogram
                    hist_matrix[i-zero_sup][j-zero_sup].Fill(pha_data['TRK_PHI2'][m]/180*TMath.Pi())
                    hist_matrix[i-zero_sup][j-zero_sup].SetTitle(name_matrix[i-zero_sup][j-zero_sup])
                    hist_matrix[i-zero_sup][j-zero_sup].GetXaxis().SetTitle("phi (rad)")
                    hist_matrix[i-zero_sup][j-zero_sup].GetYaxis().SetTitle("counts")

            hist_matrix[i-zero_sup][j-zero_sup].Draw()

            print name_matrix[i-zero_sup][j-zero_sup]

            #valore = raw_input('continue?') 

    f.Write()    

    f.Close()
