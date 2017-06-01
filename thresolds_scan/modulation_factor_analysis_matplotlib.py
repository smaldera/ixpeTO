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
import numpy as np
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

        self.c_init=0
        #self.h_modulation_factors = TH2F()

    def histogram_fitter(self):
        
        #valore = raw_input('continue?')

        #momAnalysisThresholdScan = thresholdScan() #from mom_analysis_threshold_scan_class.py
        #name_matrix = momAnalysisThresholdScan.getNameMatrix()
        name_matrix = self.name_matrix
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

        #self.h_modulation_factors = TH2F("h_modulation_factors", "h_modulation_factors", len_thresholds1, thresholds1[0]-0.5, thresholds1[len_thresholds1-1]+0.5, len_thresholds2, thresholds2[0]-0.5, thresholds2[len_thresholds2-1]+0.5)

        out_file = TFile("out_file.root", "recreate")

        #print thresholds1
        #print thresholds2

        #valore = raw_input('continue?')
        fill_matrix = []
        for i in range (0,len_thresholds1):
            fill_row = []
            #for j, t2 in enumerate(thresholds2[:6]): #j corre sull'indice, t2 sul valore
            for j, t2 in enumerate(thresholds2): #j corre sull'indice, t2 sul valore
                #DA PROVARE: usare self.hist_matrix dall'altra classe, cosi' sovrascrivo gli istogrammi e posso salvarli su ROOT
                h = inFile.Get(name_matrix[i][j])
                fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -TMath.Pi(), TMath.Pi())
                #fitFunc = TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -180, 180)  #LBaldini, presentazione per csn2, p. 11

                fitFunc.SetParLimits(0,0,100000000)     # offset >0
                fitFunc.SetParLimits(1,0,100000000)     # se c'e' bisogno di un'inversione, che la becchi con la fase [2]!!! 
                #c = TCanvas("c","c",0) #si puo' pensare di definirlo altrove
                #c.cd()
                h.Fit("fitFunc","MR")   # opzione N
                h.Draw()
                #valore = raw_input('continue?')
                out_file.cd()
                h.Write()
                #c.Close() #da togliere          #salvare gli istogrammi con fit su file root #quando fitto, bloccare per vedere un fit alla volta!!!
                modulationFactor = fitFunc.GetParameter(1)/(fitFunc.GetParameter(1)+2*fitFunc.GetParameter(0)) #NO!!! Da cambiare!!!!
                #print modulationFactor
                fill_row.append(modulationFactor)
                #self.h_modulation_factors.Fill(thresholds1[i],thresholds2[j],modulationFactor)
            fill_matrix.append(fill_row)
        #print fill_matrix
        fill_matrix = np.array(fill_matrix)
        #print fill_matrix



        fig = plt.figure(facecolor = 'white')
        ax  = fig.add_subplot(111)
        cax = ax.matshow(fill_matrix.T, origin = 'lower', aspect = 'auto', cmap = 'bone') #Spectral
        cb = plt.colorbar(cax)
        ax.set_xticklabels(['']+list(thresholds1)) #['']+list  --> con [''] mette i tick al centro dei bin
        ax.set_yticklabels(['']+list(thresholds2))
        plt.xlabel('1$^{st}$ threshold') #stringhe come LaTeX!!
        plt.ylabel('2$^{nd}$ threshold')
        #cb.set_clim(vmin=0, vmax=1)
        cb.set_label('modulation factor $\mu$', rotation = 90)
        plt.title('%s' %os.path.basename(FILE_PATH).replace('.root','')) # %.3f = float con 3 cifre dopo la virgola


        integral_thresholds1 = np.sum(fill_matrix.T, axis=1) #proiettiamo sull'asse y
        integral_thresholds2 = np.sum(fill_matrix.T, axis=0) #proiettiamo sull'asse x

        plt.figure(facecolor = 'white')
        plt.plot(thresholds2, integral_thresholds1, 'o--', label = 'thresholds 2')
        plt.plot(thresholds1, integral_thresholds2, 'o--', label = 'thresholds 1')
        plt.xlabel('threshold') #stringhe come LaTeX!!
        plt.ylabel('$\Sigma_{i}\mu_{i}$')
        plt.title('thresholds trend')
        plt.grid()
        plt.legend(fontsize = 14)

        from itertools import product
        ccc = np.argmax(fill_matrix)
        aaa = np.array(list(product(thresholds1, thresholds2)))
        print aaa[ccc]
        #print fill_matrix.T.reshape(0, len_thresholds1*len_thresholds2)[ccc]     

        #plt.savefig(....) <-- forse in conflitto con plt.show()



        '''
        self.h_modulation_factors.GetXaxis().SetTitle("1st pass threshold")
        self.h_modulation_factors.GetYaxis().SetTitle("2nd pass threshold")
        self.h_modulation_factors.Draw("colZ")
        self.h_modulation_factors.Write()
        '''

        plt.show()
    
        #valore = raw_input('continue?')
    
        out_file.Close()
        inFile.Close()
    
        #TO DO:

        #draw TH2 in 3-d??? (surf2)
        #set Stat per chi quadro (opzione di TH1/2 o canvas?? altrimenti da TBrowser View-->Editor-->Chi!!
        #istogramma del chi quadro ridotto!! 
        #che errore da' di default sull'altezza delle barre dell'istogramma
        #soglie t.c. mod_factor == max ---> draw on the TH2 --> write canvas con h_modulation_factors + TMarker
        #angolo fisso, diverse energie
        #energia fissa, diversi angoli
        #con diverso asse di polarizzazione non dovrebbe cambiare



def test():

    cc = TCanvas("cc","cc",0)

    thresholdScanObj = thresholdScan()  # oggetto della classe thresholdScan
                                        # definita in mom_analysis_threshold_scan_classe.py

    modFact = modulationFactor(thresholdScanObj)    # oggetto della classe modulationFactor
                                                    # definita in questo file
    modFact.c_init = cc
    modFact.histogram_fitter()



if __name__ == "__main__":



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
    
    print FILE_PATH
    print os.path.basename(FILE_PATH)


