#!/usr/bin/env python
# Copyright (C) 2007--2016 the X-ray Polarimetry Explorer (XPE) team.
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


#PROIEZIONE sulla DIREZIONE RICOSTRUITA

import ROOT
import numpy
from array import array
import math

import smoothing_passabassoSimo  as smooth_simo
from xpeSimo_ttree import *

from pyxpe.recon.geometry import xpePoint2d, xpeRay2d
from pyxpe.utils.logging_ import logger

import matplotlib.pyplot as plt
from pyxpe.recon.xpol import XPOL_COLUMN_PITCH
from pyxpe.recon.geometry import xpePoint2d, xpeRay2d
from pyxpe.recon.xpol import xpeHexagonCollection, adc2colors
from pyxpe.recon.clustering import hierarchical_clustering
from pyxpe.recon.recon import xpePixyRecon



class xpeSimo:

    def __init__(self, cluster, recon, raggioCut, dividiBins, baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw):

        self.event_id   = -1
        self.recon      = recon
        self.cluster    = cluster
        self.baricenter = cluster.baricenter
        self.conversion_point=recon.conversion_point
        self.phi0       = recon.phi0   			#DIREZIONE asse principale???
        self.phi1       = recon.phi1   			#DIREZIONE ricostruita???
        self.phiTang    = 0

        self.mom2_long  = recon.ma0.mom2_long		#addL
        self.mom2_trans = recon.ma0.mom2_trans		#addL
        self.mom3_long  = recon.mom3_long		#addL

        self.xnew       = 0
        self.ynew       = 0 		#nuovo punto di conversione!!
                             
        self.distConv   = 0.  		#distanza  punto di conversione rec standard, lungo la traccia
        self.x_picco    = 0.   		#

        # parametri configurabili...
        self.dividi_bins= dividiBins	# scala il n. di bins per l'istogramma: 
        self.raggioCut  = raggioCut     # taglia i pixel piu' lontani di r dalla traccia fittata ... non sembra molto utile!! 
        self.baryPadding= baryPadding   # distanza permessa tra traccia e punto di conversione analisi standard
        self.peakFinding= findMaxAlg    # seleziona l'algoritmo per la ricerca picco auger.  
        self.pcubo      = pcubo         # 0 parabola, 1 cubo
        self.draw       = draw          # disgna singoli eventi

        # parametri algoritmo peakfinding di TSpectrum
        self.maxnP=maxnP
        self.Psigma=Psigma
        self.Pthr=Pthr

        # definizione istogrammi distribuzione longitudinale carica
        self.x1=-2 
        self.x2=2 
        self.y1=-2
        self.y2=2

        # n. di bin per istogrammi!!
        self.nbinsX=int((self.x2-self.x1)/0.01);
        self.nbinsY=int((self.y2-self.y1)/0.01)

        # istogrammi
        self.h1= ROOT.TH1F("h1","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)     

        #self.h1L= ROOT.TH1F("h1L","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)
        #self.h1L_smoothSimo= ROOT.TH1F("h1L_smoothSimo","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)
        #self.h1L_ave= ROOT.TH1F("h1L_ave","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)

       
        # funzione fit istogramma long. cariche
        self.fFit_histo=0
        self.redChi2=0
        self.sumpars=array('d',[0.]*7) #????????????????????????????????????????????
        
        #variabile cui passo canvas in --> test
        self.c_init=0

        #legenda scatter plot
        self.legend = ROOT.TLegend(.1,.7,.4,.9,'   Legend') #add options

        #output ROOT file
        self.outRootFile=ROOT.TFile()


        self.bary1=ROOT.TMarker()
        self.bary2=ROOT.TMarker()
        self.line1=ROOT.TF1()
        self.line2=ROOT.TF1()
        #self.f_p3=ROOT.TF1()
        self.h2 = ROOT.TH2F()
        self.profx=ROOT.TProfile()
        self.newPoint=[0,0]
        self.newPointMarker=ROOT.TMarker() 	#addL -> tolto da funzione draw_simo e definito come oggetto della classe -> elimina problemi sul draw

        self.ellipse=ROOT.TEllipse()		#addL
        
        self.convPoint = ROOT.TMarker()		#addL per l'istogramma a dx
        self.peakPoint = ROOT.TMarker()		#addL per l'istogramma a dx


    '''
    def fit_spline(self, zero_suppression=12):
        VEDI rec_simo.py

    def fit_spine(self, num_nodes=5):
        VEDI rec_simo.py
    '''

        
    def draw(self, coordinate_system, color_map='Reds', hexcol_padding=0.1,  text=True, show=True):
        """Draw the cluster. To be moved into a separate module.
        """
        hit_positions = numpy.vstack((self.cluster.x, self.cluster.y),).transpose()
        colors = adc2colors(self.cluster.adc_values, 0, color_map)
        if coordinate_system == 'pixy':
            angle = numpy.pi/2.
        else:
            angle = 0
        hex_col = xpeHexagonCollection(padding=hexcol_padding,
                                       offsets=hit_positions, rotation=angle,
                                       edgecolors='gray', facecolors=colors)
        fig = hex_col.figure
        if text:
            adc_ref = 0.5*self.cluster.adc_values.max()
            for x, y, val in zip(self.cluster.x, self.cluster.y, self.cluster.adc_values):
                if val < adc_ref:
                    col = 'black'
                else:
                    col = 'white'
                plt.text(x, y, '%s' % val, horizontalalignment='center',
                         verticalalignment='center', size=8, color=col)
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        self.baricenter.draw(color='k')
        axis0 = xpeRay2d(self.baricenter, self.phi0)
        axis0.draw()

        
        self.conversion_point.draw(color='b')
        axis1 = xpeRay2d(self.conversion_point, self.phi1)
        axis1.draw(color='b')      
#        self.fit_spline()
#        self.fit_spine()

        if show:
            plt.show()
        return fig


    
    def rec_simo(self, coordinate_system, color_map='Reds', hexcol_padding=0.1,  text=True, show=True):

        print ""
        print "==================== "
        print ""
        
        import time
        import math
        print "SIMO: inizio draw_root"
                
        x0=self.conversion_point.x()
        y0=self.conversion_point.y()
        phi0=self.phi0  #NEW
        phi1=self.phi1  #angolo di rotazione asse ricostruito??? #NEW (era: phi = self.phi1)
        phi = phi1      #cosi' posso lasciare alcune cose scritte da simo usando phi

        print 'phi0 = ', phi0
        print 'phi1 = ', phi1

        dx = (self.cluster.x - x0)      #x0 e y0 sono le coord del punto attorno a cui ruoto (conversion point ottenuto dall'analisi standard: OK!)
        dy = (self.cluster.y - y0)
        

        #ruoto il SR di phi0 (addLori LORI): asse principale = asse x
        #ARRAY di COORDINATE
        xp = numpy.cos(phi0)*dx + numpy.sin(phi0)*dy    #x nel nuovo SR con asse ricostruito = asse x
        yp = -numpy.sin(phi0)*dx + numpy.cos(phi0)*dy   #y nel nuovo SR

        #angoli phi0 e phi1 nel nuovo SR
        phi0_bary=self.phi0-phi0      #angolo dell'asse principale (sempre 0!)
        phi1_bary=self.phi1-phi0      #angolo della direzione ricostruita
        #le loro tangenti sono i coeff. angolari delle due rette
        
        print 'phi0_bary', phi0_bary
        print 'phi1_bary', phi1_bary # == 0 perche' asse x = direz. ricostuita

        maxX= max(xp)
        minX=min(xp)
        maxY= max(yp)
        minY=min(yp)
        
        #print "DRAW SIMO: minX= ",minX, "maxX = ",maxX
        

        #Scatter plot dei pixel
        self.h2 = ROOT.TH2F("h2","",self.nbinsX,self.x1,self.x2,self.nbinsY,self.y1,self.y2)
        self.profx=ROOT.TProfile("profx","profile",100,-2,2,"")		#a cosa serve???


        # probabilmente inutile, giusto per essere sicuro che siano vuoti
        self.h1.Reset()
        #self.h1L.Reset()
        #self.h1L_smoothSimo.Reset()
        #self.h1L_ave.Reset()
        self.h2.Reset()
        self.profx.Reset()

        
        for i in range (0,len(self.cluster.x)):
             self.h2.Fill(xp[i], yp[i],self.cluster.adc_values[i])	#fill scatter plot

             self.profx.Fill(xp[i], yp[i],  numpy.sqrt(self.cluster.adc_values[i]))
									# peso con la radice dei conteggi... assumo che siano sempre >0!!!
        
             #histogram: ADC counts --> proiezione lungo asse ricostruito
             self.h1.Fill(xp[i]-minX,self.cluster.adc_values[i])	#traslazione: l'estremo sx dell'istogramma sia a zero
             self.h1.SetAxisRange(minX-minX-0.2, maxX-minX+0.2)
             self.h1.SetTitle("Proiezione lungo asse principale")

        self.h2.GetXaxis().SetRangeUser(-0.7,0.7)
        self.h2.GetYaxis().SetRangeUser(-0.6,0.6)



        #SR con ASSE PRINCIPALE dell'ellisse d'inerzia == ASSEX (fatto da LORI)
        #calcolo coordinate bary1 (baricentro) e 2 (conversion point) nel sistema roto traslato:

        #BARICENTRO del CLUSTER = bary1
        dx_bary1=(self.baricenter.x()-x0)
        dy_bary1=(self.baricenter.y()-y0)
        x_bary1= numpy.cos(phi0)*dx_bary1 + numpy.sin(phi0)*dy_bary1
        y_bary1= -numpy.sin(phi0)*dx_bary1 + numpy.cos(phi0)*dy_bary1

        self.bary1=ROOT.TMarker(x_bary1,y_bary1,20)
        self.bary1.SetMarkerColor(2) #kRed == 2


        #CONVERSION point (from standard recon)
        dx_bary2=(self.conversion_point.x()-x0) # = 0 Il punto di conversione nella ricostruzione standard e' (0,0) nel nostro SR!
        dy_bary2=(self.conversion_point.y()-y0) # = 0
        x_bary2= numpy.cos(phi0)*dx_bary2 + numpy.sin(phi0)*dy_bary2
        y_bary2 = -numpy.sin(phi0)*dx_bary2 + numpy.cos(phi0)*dy_bary2
        
        self.bary2=ROOT.TMarker(x_bary2, y_bary2,20)
        self.bary2.SetMarkerColor(4)

        #ELLIPS of INERTIA
        self.ellipse = ROOT.TEllipse(x_bary1,y_bary1,math.sqrt(self.mom2_long), math.sqrt(self.mom2_trans), 0,360, 0) #self.phi0*ROOT.TMath.RadToDeg())


        #PRINCIPAL AXIS parameters --> line1
        m=(math.tan(phi0_bary))
        q=y_bary1-m*x_bary1

        self.line1=ROOT.TF1("f1","[0]*x+[1]",minX, maxX)
        self.line1.SetParameters(m,q)
        self.line2.SetLineColor(2) #red


        #RECON DIRECTION parameters --> line2
        m2=(math.tan(phi1_bary)) #line2, blu, = 0
        q2=y_bary2-m2*x_bary2

        self.line2=ROOT.TF1("f2","[0]*x+[1]",minX,maxX)
        self.line2.SetParameters(m2,q2)
        self.line2.SetLineColor(4) #blue

        fPrincipalAxis=ROOT.TF1("fPrincipalAxis"," [0] + [1]*x", minX, maxX)
        fPrincipalAxis.SetParameter(0,q)
        fPrincipalAxis.SetParameter(1,m)
        
        par = [q2,m2]
        
        '''
        self.f_p3=ROOT.TF1("f_p3","[0]*x*x*x+[1]*x*x+[2]*x+ ([3]- [0]*[4]*[4]*[4]-[1]*[4]*[4]-[2]*[4] )", minX, maxX) # pol3  per bary2 
        # fissa il passaggio per il punto di conv. della rec standard
        #self.f_p3.FixParameter(3,y_bary2)
        #self.f_p3.FixParameter(4,x_bary2)
        .
        .
        (vedi rec_simo.py)
        .
        .
        #self.h1L_smoothSimo=self.smooth_simo(self.h1L)
        #self.h1L_ave=self.media_mobile(self.h1L)
        '''

        #SMOOTHING histogram h1 (ADC counts --> proiezione lungo asse ricostruito)
        self.h1_smoothSimo=self.smooth_simo(self.h1)

      
        #GET Auger Peak point from histogram and get its coordinates on scatter plot
        self.newPoint= self.cerca_piccoAugerElectron(self.h1_smoothSimo, fPrincipalAxis, minX, par) #addL: need to project on principal axis
											
        print "newPoint (PeakImpP) x= ", self.newPoint[0], "y = ", self.newPoint[1]  	#addL
        self.newPointMarker=ROOT.TMarker(self.newPoint[0], self.newPoint[1],34) 	#addL
        self.newPointMarker.SetMarkerColor(2) #kRed == 2                        	#addL
        self.newPointMarker.SetMarkerSize(1)                                    	#addL
        #valore = raw_input('aaaaaa????')
              
        new_coord=self.undo_rotoTraslation (self.newPoint[0],self.newPoint[1])
        self.xnew=new_coord[0]
        self.ynew=new_coord[1]
        print "AAA x= ", self.xnew, "y = ", self.ynew
        




    # DRAW function
    def draw_simo(self):

        #CANVAS
        self.c_init.Clear()		#canvas passed to c_init in --> test
        self.c_init.Divide(2,1)

	#SX SUB CANVAS
        self.c_init.cd(1)
 
	#PIXELS SCATTER PLOT
        self.h2.Draw("colZ")

        self.profx.Draw("samep")	#????????????????????

	#ELLIPSE of INERTIA
        #ellipse=ROOT.TEllipse(self.baricenter.x(),self.baricenter.y(), math.sqrt(self.mom2_long), math.sqrt(self.mom2_trans), 0,360, self.phi0*ROOT.TMath.RadToDeg()) --> NO
        self.ellipse.SetFillStyle(0)	#addL
        self.ellipse.Draw()		#addL
        

	#OTHER PLOTS
        self.bary1.Draw()		#draw barycenter
        self.bary2.Draw()		#draw conversion point
           
        self.newPointMarker.Draw()      #draw Auger electron peak along recon direction	#addL

        self.line1.Draw("samel")	#draw principal axis
        self.line2.Draw("samel")	#draw recon direction

        #self.f_p3.Draw("samel")


        #LEGEND
        self.legend.AddEntry(self.bary1, "barycenter", "p")
        self.legend.AddEntry(self.bary2, "conversion point", "p")
        self.legend.AddEntry(self.newPointMarker, "Auger peak", "p")
        self.legend.AddEntry(self.line1, "principal axis", "l")
        self.legend.AddEntry(self.line2, "recon direction", "l")
        self.legend.Draw("same")


        #DX SUB CANVAS
        self.c_init.cd(2)
        self.h1.Draw("hist")
	
        #SMOOTHED HISTOGRAM
        self.h1_smoothSimo.SetLineColor(4)
        self.h1_smoothSimo.SetLineWidth(2)
        self.h1_smoothSimo.Draw("sames")
        #valore = raw_input('continue?') 

	#CONVERSION POINT
        self.convPoint = ROOT.TMarker(self.distConv,0,22)
        self.convPoint.SetMarkerColor(2)
        self.convPoint.Draw()

	#AUGER PEAK
        self.peakPoint = ROOT.TMarker(self.x_picco,0,22)
        self.peakPoint.SetMarkerColor(4)
        self.peakPoint.Draw()  
            
        #self.h1L_ave.SetLineColor(8)
        #self.h1L_ave.SetLineWidth(4)
        #self.h1L_ave.Draw("sames")


        if self.peakFinding==3:                         #???????????????????
            self.fFit_histo.Draw("samel")               #???????????????????

        if self.peakFinding==2:                         #???????????????????
            self.fFit_histo.Draw("samel")               #???????????????????
            
        self.c_init.Update() #!!!!!!!!!!!!!!!!!!!!!!!!!
        self.outRootFile.cd()
        self.c_init.Write()

        #valore = raw_input('continue?')
        
        return 1



    def deleteHistos(self):
        del  self.h1
        del  self.h1_smoothSimo     
        #del  self.h1L
        #del  self.h1L_smoothSimo
        #del  self.h1L_ave
        del  self.h2
        del  self.profx


    '''
    def min_dist (self, pixel_x, pixel_y, par, xmin,xmax): 
        """
        #cerca la minima distanza tra il punto e la funzione...
        #restituisce la x del punto sulla funzione a cui ho il minimo
        """
        (vedi rec_simo.py)
    '''      
    
    '''    
    def get_xDist (self,fLin, x_min, l):
        """ cerco l'estremo superiore tale che l'integrale tra xmin e b sia ==l   
             SAREBBE DA CALCOLARE ANALITICAMENTE!!!!!!!, ma per ora mi accontento di calcolo numerico... 
        """
        (vedi rec_simo.py)
    '''
    

    
    def smooth_simo (self, h1):
       # ficco i valore dellh su un vettore
       nbins=h1.GetNbinsX()
       x=numpy.array([0.]*nbins, dtype=float)
       for i in range (1,nbins):
           x[i]=h1.GetBinContent(i)

       print "creo oggetto..."
       simoFilter=smooth_simo.filtra_segnale()

       print "init filtro.. "
       aa=simoFilter.init_filtro();
       y=simoFilter.filtra(x)
       hF=h1.Clone()
       hF.Reset()
       for i in range (1,nbins):
           hF.SetBinContent(i,y[i])
           #print "filtered: i= ",i," y = ",y[i]

       return hF



    def media_mobile (self, h1):        
       nbins=h1.GetNbinsX()
       hF=h1.Clone()
       hF.Reset()
       for i in range (2,nbins-2):
             
             ave=  (h1.GetBinContent(i)+ h1.GetBinContent(i-1)+h1.GetBinContent(i+1))/3.
             hF.SetBinContent(i,ave)
       return hF      



    def undo_rotoTraslation (self,xp,yp): #modified by Lori: here phi0, not phi1 (principal axis, not recon direction)

       #phi=self.phi0
       #x0=self.baricenter.x()
       #y0=self.baricenter.y()
       phi0=self.phi0
       x0=self.conversion_point.x()
       y0=self.conversion_point.y()  

       dx=xp*numpy.sin(phi0)+yp*numpy.cos(phi0)
       dy=xp*numpy.cos(phi0)-yp*numpy.sin(phi0)

       x=dx+x0  
       y=dy+y0
       new_coord=[x,y]
       return new_coord        

   

    def fit2Gaussiane (self, h1):

       G0 = ROOT.TF1 ("G0","gaus",0.2,0.6)
       h1.Fit("G0","MER")
       mean0=G0.GetParameter(1)
       sigma0=G0.GetParameter(2)

       G2 = ROOT.TF1 ("G2","gaus",0.25,0.5)
       G2.SetParameter(0, 9.74680e+02 )
       G2.SetParameter(1,  3.17394e-01)
       G2.SetParameter(2, 5.47369e-02 )

       G2.SetParLimits(1, mean0-sigma0/2., mean0+sigma0/2.)
       #G2.FixParameter(1)
       #G2.FixParameter(2)
       
       h1.Fit("G2","MER")
        
       G1 = ROOT.TF1 ("G1","gaus",0.05,0.35)
       G1.SetParameter(0,3.33294e+02 )
       G1.SetParameter(1, 1.35906e-01 )
       G1.SetParameter(2, 5.30043e-02  )

       G1.SetParLimits(1, 0.03,0.25)
       G1.SetParLimits(2, 0.02,0.1)
       
       #G1.FixParameter(1)
       #G1.FixParameter(2)
       
       h1.Fit("G1","MER")


       Gsum=ROOT.TF1("Gsum","gaus(0)+gaus(3)",-0.05,0.7)
       Gsum.SetParameter(0, G1.GetParameter(0))
       Gsum.SetParameter(1, G1.GetParameter(1))
       Gsum.SetParameter(2, G1.GetParameter(2))
       Gsum.SetParameter(3, G2.GetParameter(0))
       Gsum.SetParameter(4, G2.GetParameter(1))
       Gsum.SetParameter(5, G2.GetParameter(2))

       Gsum.SetParLimits(1, 0.03,0.3)
       Gsum.SetParLimits(2, 0.02,0.1)
       Gsum.SetParLimits(4,  mean0-sigma0/2., mean0+sigma0/2.)
       
            
       Gsum.SetLineColor(2)
       h1.Fit("Gsum","MER")
               
       return Gsum
       


    def fitGaus_cutoff (self, h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print "fit con cut-off"
        
        max = h.GetBinCenter( h.GetMaximumBin())
        convp=self.distConv
        print "distConv= ",self.distConv

         
        self.c_init.cd(2)
        #cccc = ROOT.TCanvas("cccc", "cccc", 0) #addL
        #cccc.cd()        #addL
        G0 = ROOT.TF1("G0","gaus",max-0.1,max+0.1)
        h.Fit("G0","LMER")
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        G0.Draw("samel")
        
         
        #OK !!
        #G1=ROOT.TF1 ("G1","gaus",0.05,mean0-sigma0)
        G1=ROOT.TF1 ("G1","gaus",convp-0.1,convp+0.1)
        
        G1.SetLineColor(4)
        h.Fit("G1","LMER")
        G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)
                
        #ok2
        #f1= ROOT.TF1 ("f1","  ([0]*x-[1])*(1./  (1+exp( (x-[2])/[3])  ))",mean1+sigma1,max+0.2)
        #f1.SetParameters(500,88,mean0,0.01)
        f1= ROOT.TF1 ("f1","  ([0]*x-[1])*(1./  (1+exp( (x-[2])/[3])  ))", mean0-sigma0,max+0.2)
        f1.SetParameters(500,88,mean0,0.01)

        
        f1.SetParLimits(2, mean0-sigma0/2., mean0+sigma0)
        #f1.FixParameter(2, mean0)
        h.Fit("f1","LMER")
        f1.SetLineColor(3)
        f1.Draw("samel")

        
        #OK
        fsum=ROOT.TF1("fsum","gaus(0)+ (  ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))*(  ( ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))) >0  )  )",convp-0.1,max+0.2)

        #fsum=ROOT.TF1("fsum","landau(0)+ (  ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))*(  ( ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))) >0  )  )",convp-0.1,max+0.2)

                
        #ok
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        
        fsum.SetParameter(0, 200)
        #fsum.SetParameter(1, max-sigma0 )
        #fsum.SetParameter(2, sigma0 ) #!!!!!!!!!!!!!!!

        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)

        fsum.SetParLimits(0, 50,300)
        #fsum.SetParLimits(1, 0.02,max-sigma0)
        #fsum.SetParLimits(1, convp-0.1 , min(convp+0.1,max) ) 
        fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0/2.) ) 
        
        fsum.SetParLimits(2, 0.02, sigma0)
        
        #fsum.SetParLimits(5,  mean0-2.*sigma0, mean0+2.*sigma0)
        #fsum.FixParameter(5,  0.237)
        
        #OK
        fsum.SetParLimits(5,  mean0-sigma0, mean0+sigma0)

        fsum.SetLineColor(2)
        fsum.SetLineWidth(4)
        h.Fit("fsum","LMER")

        self.redChi2=fsum.GetChisquare()/fsum.GetNDF()
        #fsum.Draw("samel")

        #cccc.Close()

        return fsum
   



    #per SR con asse X = asse ricostruito
    #def cerca_piccoAugerElectron (self, h, fLin, minX,par):

    def cerca_piccoAugerElectron(self, h, fPrincipalAxis, minX, par): #modified by LORI

    #def cerca_piccoAugerElectron(self, h, minX):
       # h e' l'istogramma in funzione della lunghezza lungo la traccia
       # restituisce le coord del picco auger 

        n_foundPeaks=0

        x_picco=0
      
        print 'peakFinding = ', self.peakFinding
        #valore = raw_input('continue?')  

        #cerca picchi!
        if (self.peakFinding==1):
            print " inizio ROOT.Tspectrum"
            s=ROOT.TSpectrum(self.maxnP)
            n_foundPeaks=s.Search(h,self.Psigma,"",self.Pthr)               
            x_peaks_raw=s.GetPositionX()
            x_peaks=[]
        
            for i in range (0,n_foundPeaks):
                x_peaks.append(x_peaks_raw[i])

            x_peaks.sort()
 
            #for i in range (0,n_foundPeaks):
            #   print "x peak found = ",x_peaks[i] 
        
            if n_foundPeaks >=2:
                self.x_picco=x_peaks[0]
            else:
                newPoint=[-100,100]
                return newPoint
        #------------------------        
                
        if (self.peakFinding==2): 
            # uso fit doppia gaussiana
            self.fFit_histo=self.fitGaus_cutoff(h) 
            self.x_picco=self.fFit_histo.GetParameter(1)


        if (self.peakFinding==3):
            self.fFit_histo=self.fit2Gaussiane(h)
            self.x_picco=self.fFit_histo.GetParameter(1)
            
         # converto allo spazio xy    
                    
        x_lin= minX + self.x_picco #-self.conversion_point.x()
        #par = [q,m]
        y_lin= par[0] + par[1]*x_lin  #addL   #?????????????????????????????????????????????????????????????????????????????????????????
			#che y attribuisco? quella della retta della direzione ricostruita? o quella di una spline valutata in quel punto?
                       	#per ora gli ho dato la y presa dalla retta della direzione ricostruita, che qui NON e' l'asse x.	

          
        print 'minX = ', minX
        print 'x_Picco (hist) = ', self.x_picco 
        print 'x_lin = ', x_lin
        print 'y_lin = ', y_lin

        newPoint=[x_lin,y_lin]
        return newPoint
         

  

def test(filePath, num_events,raggioCut, dividiBins, baryPadding, findMaxAlg,  zero_suppression=9, coordinate_system='xpedaq', pcubo=0, maxnP=4, Psigma=2, Pthr=0.0001, draw=False):
       """
       """

       print "raggioCut = ",raggioCut
       print "dividiBins= ",dividiBins
       print "baryPadding = ",baryPadding 
       print "findMaxAlg = ",findMaxAlg
       print "pcubo = ", pcubo
       print "maxnP = ", maxnP
       print "Psigma ",Psigma
       print "Pthr = ",Pthr

       
       h_phi1=ROOT.TH1F("h_phi1","",360,0,360)
       h_phi_tang=ROOT.TH1F("h_phi1_tang","",360,0,360)


       h_x=ROOT.TH1F("h_x","",200,-0.5,0.5)
       h_x1=ROOT.TH1F("h_x1","",200,-0.5,0.5)

       h_y=ROOT.TH1F("h_y","",200,-0.5,0.5)
       h_y1=ROOT.TH1F("h_y1","",200,-0.5,0.5)



       #addL
       c= ROOT.TCanvas("c","c", 2000,1000)
       c.Divide(2,1) 


       #if salva tree?
       event_id=0
       outRootFile=ROOT.TFile("out.root","recreate")
       myTree=myTTree()  # from xpeSimo_ttree.py
                           
       
       
       from pyxpe.recon.binio import xpeBinaryFileWindowed
       input_file = xpeBinaryFileWindowed(filePath)
       for i in xrange(num_events):
           event = input_file.next()
           #print event
           #print "event id = ",event_id
           cluster = hierarchical_clustering(event, zero_suppression, coordinate_system)[0]

           
           print cluster
           recon = xpePixyRecon(cluster)
           if not recon.error_summary:             
               xpeSimoAA=xpeSimo(cluster,recon,raggioCut,dividiBins,baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw)
               xpeSimoAA.c_init=c    #addL. passo un canvas per poter disegnare sempre sullo stesso (sicuramente c'e' un modo piu' furbo!!!!! )
               xpeSimoAA.event_id=event_id
               xpeSimoAA.outRootFile= outRootFile
                   
                              
               recSimo=xpeSimoAA.rec_simo(coordinate_system)
               if draw:
                   xpeSimoAA.draw_simo()
                   valore = raw_input('continue?') 
               myTree.Fill(xpeSimoAA)
               
                                           
               
               h_phi1.Fill(xpeSimoAA.phi1*ROOT.TMath.RadToDeg() )
               h_phi_tang.Fill(xpeSimoAA.phiTang*ROOT.TMath.RadToDeg() )
               if (xpeSimoAA.xnew!=-100):
                   h_x.Fill(xpeSimoAA.conversion_point.x())
                   h_x1.Fill(xpeSimoAA.xnew)
                   h_y.Fill(xpeSimoAA.conversion_point.y())
                   h_y1.Fill(xpeSimoAA.ynew)

               xpeSimoAA.deleteHistos()
           event_id =event_id+1
           #valore = raw_input('continue?')               
     
       c3=ROOT.TCanvas("c3","",2000,1000)
       c3.Divide(1,2)
       

       #h_phi1.Draw()
       #h_phi_tang.SetLineColor(2)
       #h_phi_tang.Draw("sames")
       h_x1.SetLineColor(2)
       h_y1.SetLineColor(2) 
       c3.cd(1)
       h_x.Draw()
       h_x1.Draw("sames")
       c3.cd(2)
       h_y.Draw()
       h_y1.Draw("sames")


       
       # scrivi outfile 
       miofile = open('miofile.txt','w')   
       miofile.write(str(raggioCut)+ " "+str(dividiBins)+ "  "+str(baryPadding)+"  "+str(findMaxAlg)+" "+str(pcubo)+ "  "+str(maxnP)+"  "+str(Psigma)+" "+str(Pthr)+"  " +str( h_x1.GetRMS() )+" \n")

       #outRootFile=ROOT.TFile("out.root","recreate")
       outRootFile.cd()
       h_phi1.Write()
       h_phi_tang.Write()
       h_y.Write()
       h_y1.Write()
       h_x.Write()
       h_x1.Write()
       c3.Write()

       myTree.treeSimo.Write()
       
       outRootFile.Close()
       

       
       valore = raw_input('continue?')    
        
if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-n', '--num_events', type=int, default=100,
                        help = 'number of events to be processed')
    parser.add_argument('-z', '--zero-suppression', type=int, default=9,
                        help = 'zero-suppression threshold')
    parser.add_argument('-c', '--coordinate-system', type=str, default='pixy',
                        help = 'coordinate system for the clustering')


    parser.add_argument('-r', '--raggioCut', type=float, default=0.07,
                        help = 'raggio intorno al fit per accetare i pixel da proiettare')


    parser.add_argument('-divbins', '--dividiBins', type=float, default=0.5,
                        help = 'fattore per n. di bins histo proiettato')


    parser.add_argument('-baryPadding', '--baryPadding', type=float, default=0.005,
                        help = 'limite distanxa funzione da baricentro ')


    parser.add_argument('-findMaxAlg', '--findMaxAlg', type=int, default=1,
                        help = 'algoritmo ricerca picco e auger... 1 TSpectrum, 2->due gauss - 3-> fit gaus + cutoff ') 
     

    parser.add_argument('-pcubo', '--pcubo', type=int, default=0,
                        help = 'se 0 fissa il parametro di 3 grado a zero -> parabola!  ')
    

    parser.add_argument('-d', '--draw', type=bool, default=False,
                        help = 'draw (da aggiungere storage  su file)  (True/ False)  ')
    
    parser.add_argument('-maxnP', '--maxnP', type=int, default=4,
                        help = 'max numner of peaks in TSectrum constructor  ')

    parser.add_argument('-Psigma', '--Psigma', type=int, default=2,
                        help = 'sigma in TSectrum search peak   ')

    parser.add_argument('-Pthr', '--Pthr', type=float, default=0.0001,
                        help = 'threshold in TSectrum search peak   ')
    
    args = parser.parse_args()
    test(args.infile, args.num_events, args.raggioCut,  args.dividiBins, args.baryPadding, args.findMaxAlg, args.zero_suppression,
         args.coordinate_system, args.pcubo, args.maxnP, args.Psigma, args.Pthr,  args.draw )
