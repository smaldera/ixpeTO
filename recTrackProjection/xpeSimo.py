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


import ROOT
import numpy
from array import array
import math


import smoothing_passabassoSimo  as smooth_simo
from xpeSimo_ttree import *


import matplotlib.pyplot as plt
from gpdswswig.Recon import *

class xpeSimo:

    def __init__(self, track,   raggioCut,dividiBins, baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw  ):

        self.event_id=-1
        self.baricenter_X=0
        self.baricenter_Y=0
        self.conversion_point_X=0
        self.conversion_point_Y=0
       
        self.phi0=0
        self.phi1=0
        self.phiTang=0

        self.xnew=0
        self.ynew=0 #  nuovo punto di conversione!!
               
              
        self.distConv=0.  #distanza  punto di conversione rec standard, lungo la traccia
        self.x_picco=0.   #

        # parametri configurabili...
        self.dividi_bins=dividiBins       # scala il n. di bins per l'istogramma: 
        self.raggioCut=raggioCut          # taglia i pixel piu' lontani di r dalla traccia fittata ... non sembra molto utile!! 
        self.baryPadding=baryPadding      # distanza permessa tra traccia e punto di conversione analisi standard
        self.peakFinding=findMaxAlg       # seleziona l'algoritmo per la ricerca picco auger.  
        self.pcubo=pcubo                  # 0 parabola, 1 cubo
        self.draw=draw                    # disgna singoli eventi

        # parametri algoritmo peakfinding di TSpectrum
        self.maxnP=maxnP
        self.Psigma=Psigma
        self.Pthr=Pthr

        # definizione  istogrammi distribuzione longitudinale carica
        self.x1=-2 
        self.x2=2 
        self.y1=-2
        self.y2=2

        # n. di bin per istogrammi!!
        self.nbinsX=int( (self.x2-self.x1)/0.01);
        self.nbinsY=int((self.y2-self.y1)/0.01)

        # istogrammi
        self.h1= ROOT.TH1F("h1","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)     
        self.h1L= ROOT.TH1F("h1L","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)
        self.h1L_smoothSimo= ROOT.TH1F("h1L_smoothSimo","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)
        self.h1L_ave= ROOT.TH1F("h1L_ave","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)

        
      

        

        # funzione fit istogramma long. cariche
        self.fFit_histo=0
        self.redChi2=0
        self.sumpars=array('d',[0.]*7) 
        
        self.c_init=0
        # ROOT.TCanvas("cc","cc", 2000,1000) 
        #self.c_init.Divide(2,1)
        #self.c_init.Draw()
        self.outRootFile=ROOT.TFile()

        self.bary1=ROOT.TMarker()
        self.bary2=ROOT.TMarker()
        self.line1=ROOT.TF1()
        self.line2=ROOT.TF1()
        self.f_p3=ROOT.TF1()
        self.h2 = ROOT.TH2F()
        self.profx=ROOT.TProfile()
        self.newPoint=[0,0]
        
       
 """
    def fit_spline(self, zero_suppression=12):
        """To be moved into recon.
        """
        from scipy.interpolate import UnivariateSpline
        _mask = self.cluster.adc_values >= zero_suppression
        x = self.cluster.x[_mask]
        y = self.cluster.y[_mask]
        adc_values = self.cluster.adc_values[_mask]
        weights = (adc_values/float(adc_values.max()))**0.5
        dx = (x - self.baricenter.x())
        dy = (y - self.baricenter.y())
        xp = numpy.cos(self.phi0)*dx + numpy.sin(self.phi0)*dy
        yp = -numpy.sin(self.phi0)*dx + numpy.cos(self.phi0)*dy
        #s = UnivariateSpline(xp, yp, w=weights, s=0.5)
        s = UnivariateSpline(xp, yp, w=weights, s=1)
       
        _xp = numpy.linspace(xp.min(), xp.max(), 25)
        _yp = s(_xp)
        dx = numpy.cos(-self.phi0)*_xp + numpy.sin(-self.phi0)*_yp
        dy = -numpy.sin(-self.phi0)*_xp + numpy.cos(-self.phi0)*_yp
        x = dx + self.baricenter.x()
        y = dy + self.baricenter.y()
        plt.plot(x, y, '-', lw=2, color='black')
"""
        
       
    def rec_simo(self):

        print ""
        print "==================== "
        print ""
        
        import time
        import math
        #import ROOT

        #creo liste con x,y,z di ogni pixel!!! forse meglio spostarlo fuori da questa classe e passare le liste!
        n_hits=self.track.numHits()
        hit=track.hits()
        x=[0].*n_hits
        y=[0].*n_hits
        adc=[0].*n_hits
        
        for i in range (0,n_hits):
            x[i]=hit[i].x
            y[i]=hit[i].y
            adc[i]=hit[i].pulseHeight
            
        
        # old: self.baricenter, phi0
        
        x0=self.conversion_point_X
        y0=self.conversion_point_Y
        phi=self.phi0
        
        dx = (x - x0)
        dy = (y - y0)
        xp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
        yp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy

        phi0_bary=self.phi0-phi
        phi1_bary=self.phi1 -phi
               
                
        maxX= max(xp)
        minX=min(xp)
        maxY= max(yp)
        minY=min(yp)
        
        #print "DRAW SIMO: minX= ",minX, "maxX = ",maxX
        
        self.h2 = ROOT.TH2F("h2","",self.nbinsX,self.x1,self.x2,self.nbinsY,self.y1,self.y2)
        self.profx=ROOT.TProfile("profx","profile",100,-2,2,"")

        # probabilmente inutile, gisuto per essere sicuro che siano vuoti
        self.h1.Reset()
        self.h1L.Reset()
        self.h1L_smoothSimo.Reset()
        self.h1L_ave.Reset()
        self.h2.Reset()
        self.profx.Reset()
       

        
        for i in range (0,len(x)):
             self.h2.Fill(xp[i], yp[i],adc[i])
             self.profx.Fill(xp[i], yp[i],  numpy.sqrt(adc[i])) # peso con la radice dei conteggi... assumo che siano sempre >0!!!
             #h1.Fill(xp[i]-minX,self.cluster.adc_values[i])
             self.h1.Fill(xp[i]-minX,adc[i])
        
        self.h2.GetXaxis().SetRangeUser(-0.8,1)
        self.h2.GetYaxis().SetRangeUser(-0.6, 0.6)

               
        #gr1 = ROOT.TGraph (len(xp),xp,yp )

        #calcolo coodrdinate bary1 e 2 nel sistema roto traslato:

        dx_bary1=(self.baricenter_X-x0)
        dy_bary1=(self.baricenter_Y-y0)
        x_bary1= numpy.cos(phi)*dx_bary1 + numpy.sin(phi)*dy_bary1
        y_bary1 = -numpy.sin(phi)*dx_bary1 + numpy.cos(phi)*dy_bary1


        self.bary1=ROOT.TMarker(x_bary1,y_bary1,20)
        self.bary1.SetMarkerColor(2)

        dx_bary2=(self.conversion_point_X-x0)
        dy_bary2=(self.conversion_point_Y-y0)
        x_bary2= numpy.cos(phi)*dx_bary2 + numpy.sin(phi)*dy_bary2
        y_bary2 = -numpy.sin(phi)*dx_bary2 + numpy.cos(phi)*dy_bary2
        
        self.bary2=ROOT.TMarker(x_bary2, y_bary2,20)
        self.bary2.SetMarkerColor(4)

                
        m=(math.tan(phi0_bary))
        q=y_bary1-m*x_bary1

        m2=(math.tan(phi1_bary))
        q2=y_bary2-m2*x_bary2
      
        
        self.line1=ROOT.TF1("f1","[0]*x+[1]",minX, maxX)
        self.line1.SetParameters(m,q)

        self.line2=ROOT.TF1("f2","[0]*x+[1]",minX,maxX)
        self.line2.SetParameters(m2,q2)
        self.line2.SetLineColor(4)

        
        self.f_p3=ROOT.TF1("f_p3","[0]*x*x*x+[1]*x*x+[2]*x+ ([3]- [0]*[4]*[4]*[4]-[1]*[4]*[4]-[2]*[4] )", minX, maxX) # pol3  per bary2 
        
        # fissa il passaggio per il punto di conv. della rec standard
        #self.f_p3.FixParameter(3,y_bary2)
        #self.f_p3.FixParameter(4,x_bary2)

        self.f_p3.SetParLimits(3,y_bary2-self.baryPadding ,y_bary2+self.baryPadding)
        self.f_p3.SetParLimits(4,x_bary2-self.baryPadding,x_bary2+self.baryPadding)

        if self.pcubo==0:
            self.f_p3.FixParameter(0,0) # il cubo diventa una parabola!
            #self.f_p3.FixParameter(1,0) # il cubo diventa una retta!!!!!
        print "fit con  cubo... "

        
        
        #gr1.Fit("f_p3","MER")
        self.profx.Fit("f_p3","MERw")

        par3=self.f_p3.GetParameter(3)
        par4=self.f_p3.GetParameter(4)
        
        a=self.f_p3.GetParameter(0)
        b=self.f_p3.GetParameter(1)
        c= self.f_p3.GetParameter(2)
        d= par3-a*par4**3-b*par3**2-c*par4    
        
        m=3.*a*x_bary2**3+2.*b*x_bary2+c
        self.phiTang=numpy.arctan(m)-phi
        par=[a,b,c,d] 
        
        
        x_min_dist=0
        
        fLin=ROOT.TF1("fLin","pow(  (1+ pow(3*[0]*x*x+2*[1]*x+[2],2) ) ,0.5)", minX, maxX)
        fLin.SetParameter(0,a)
        fLin.SetParameter(1,b)
        fLin.SetParameter(2,c)

        for i in range (0,len(x)):
             
             x_min_dist=self.min_dist(xp[i], yp[i], par, minX, maxX )
             y_min_dist=par[0]*x_min_dist**3+par[1]*x_min_dist**2+par[2]*x_min_dist+par[3]
             
             distL=fLin.Integral(minX,x_min_dist)          
             #print "dist x= ",x_min_dist-minX, "distL = ",distL
             radius=math.sqrt( (x_min_dist-xp[i])**2 + (y_min_dist-yp[i])**2)
             
             if radius <self.raggioCut:
                 self.h1L.Fill(distL,adc[i])
            
                 
        self.distConv=fLin.Integral(minX,0.) # in questo sitema di rif. il punto di conv e' in (0,0)         
        print "DIST conv = ",  self.distConv
                
        self.h1L_smoothSimo=self.smooth_simo(self.h1L)
        self.h1L_ave=self.media_mobile(self.h1L)
        
        
        # restituisce punto sulla traccia
        self.newPoint= self.cerca_piccoAugerElectron(self.h1L_smoothSimo, fLin, minX,par)
        m=3.*a*self.newPoint[0]**3+2.*b*self.newPoint[0]+c
        self.phiTang=numpy.arctan(m)-phi       
                
        
        new_coord=self.undo_rotoTraslation (self.newPoint[0],self.newPoint[1])
        self.xnew=new_coord[0]
        self.ynew=new_coord[1]
        print "AAA x= ", self.xnew, "y = ", self.ynew
                        


        # DRAW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def draw_simo(self):
        #if self.draw==True:

          
        self.c_init.Clear()
        self.c_init.Divide(2,1)
        self.c_init.cd(1)
        
        self.h2.Draw("colZ")
        #gr1.Draw("p")
        self.profx.Draw("samep")
        #ellips=ROOT.TEllipse(self.baricenter.x(),self.baricenter.y(), math.sqrt(self.mom2_long), math.sqrt(self.mom2_trans), 0,360, self.phi0*ROOT.TMath.RadToDeg())
        #ellips.SetFillStyle(0)
        #ellips.Draw()
        
        self.bary1.Draw()   
        self.bary2.Draw()
           
        self.line1.Draw("samel")
        self.line2.Draw("samel")

           
        self.f_p3.Draw("samel")
         
        PeakImpP=ROOT.TMarker( self.newPoint[0],  self.newPoint[1],34)
        PeakImpP.SetMarkerColor(2)    
        PeakImpP.Draw()   
                              
        self.h1L.GetXaxis().SetRangeUser(-0.4,1.5)
        self.h1L.SetLineColor(1)
        self.h1L.SetLineWidth(1)
            
        
        
        self.c_init.cd(2)
        self.h1L.Draw("hist")
        convPoint = ROOT.TMarker(self.distConv,0,22)
        convPoint.SetMarkerColor(2)
        convPoint.Draw()

        peakPoint = ROOT.TMarker(self.x_picco,0,22)
        peakPoint.SetMarkerColor(4)
        peakPoint.Draw()
            
        
        
        self.h1L_smoothSimo.SetLineColor(4)
        self.h1L_smoothSimo.SetLineWidth(2)
        self.h1L_smoothSimo.Draw("sames")
            
            
        self.h1L_ave.SetLineColor(8)
        self.h1L_ave.SetLineWidth(4)
        #self.h1L_ave.Draw("sames")

        if self.peakFinding==3:
            self.fFit_histo.Draw("samel")

        if self.peakFinding==2:
            self.fFit_histo.Draw("samel")   
               
        self.c_init.Update() #!!!!!!!!!!!!!!!!!!!!!!!!!
        self.outRootFile.cd()
        self.c_init.Write()
            
        valore = raw_input('continue?')
                
        
        return 1


    def deleteHistos(self):
        del  self.h1      
        del  self.h1L
        del  self.h1L_smoothSimo
        del  self.h1L_ave
        del  self.h2
        del  self.profx

    def min_dist (self, pixel_x, pixel_y, par, xmin,xmax): 

        """
        #erca la minima distanza tra il punto e la funzione...
        #restituisce la x del punto sulla funzione a cui ho il minimo
       """
        
        x_start=xmin
        x_stop=xmax
        n_steps=30.
        x_min=0
        min=1e8
        
        for n_iters in range(0,2):     
            min=1e8
        
            dx1=(x_stop-x_start)/n_steps
            for i in range (0, int(n_steps)):
                x=x_start+i*dx1
                y=par[0]*x*x*x+par[1]*x*x+par[2]*x+par[3]
                dist=math.sqrt( (x-pixel_x)*(x-pixel_x)+(y-pixel_y)*(y-pixel_y))
                if (dist<=min): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aggiunto <= !!!
                    min=dist
                    x_min=x

            x_start=x_min-dx1
            x_stop=x_min+dx1
        #end for n_iters

        #print "min dist = ",min, "xmin = ",x_min
        return x_min
      
        
    def get_xDist (self,fLin, x_min, l):
        """ cerco l'estremo superiore tale che l'integrale tra xmin e b sia ==l   
             SAREBBE DA CALCOLARE ANALITICAMENTE!!!!!!!, ma per ora mi accontento di calcolo numerico... 
        """
        x_max=l
        n_steps=15.
        dx=(x_max-x_min)/n_steps

        sup_min=x_min
        xfound=0
        
        for n_iter in range (0,2):

            min_d=1000
            for i in range (0, int( n_steps)):
                sup=sup_min+i*dx 
                d_calc=  math.fabs(fLin.Integral(x_min,sup)-l)
            
                if (d_calc<min_d):
                    min_d=d_calc
                    xfound=sup
            
            # setta i parametri per la prossima iterazione...           
            sup_min=xfound-dx
            sup_max=xfound+dx
            dx=(sup_max-sup_min)/n_steps        
            
            
        print "returning xfound = ",xfound
        return xfound 

    

    
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


    def undo_rotoTraslation (self,xp,yp):

       phi=self.phi0
       #x0=self.baricenter.x()
       #y0=self.baricenter.y()
       phi=self.phi1
       
       x0=self.conversion_point.x()
       y0=self.conversion_point.y()

       

       dx=xp*numpy.sin(phi)+yp*numpy.cos(phi)
       dy=xp*numpy.cos(phi)-yp*numpy.sin(phi)

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

        
        G0 = ROOT.TF1 ("G0","gaus",max-0.1,max+0.1)
        G0.SetParameter(1,max)
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
        #f1.Draw("samel")

        
        #OK
        fsum=ROOT.TF1("fsum","gaus(0)+ (  ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))*(  ( ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))) >0  )  )",convp-0.2,max+0.2)

        #fsum=ROOT.TF1("fsum","landau(0)+ (  ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))*(  ( ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))) >0  )  )",convp-0.1,max+0.2)
     
        
        #ok
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        
        
        #fsum.SetParameter(0, 200)
        #fsum.SetParameter(1, max-sigma0 )
        #fsum.SetParameter(2, sigma0 ) #!!!!!!!!!!!!!!!

        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)

       
        fsum.SetParLimits(0, 30,300)
        
        
        #fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0/2.) ) 
        fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0) ) 
        
        fsum.SetParLimits(2, 0.02, sigma0)
        #fsum.SetParLimits(2, sigma1-0.01, sigma1+0.01)
        
        #fsum.SetParLimits(5,  mean0-2.*sigma0, mean0+2.*sigma0)
        #fsum.FixParameter(5,  0.237)
        
        #OK
        fsum.SetParLimits(5,  mean0-sigma0, mean0+sigma0)

        fsum.SetLineColor(2)
        fsum.SetLineWidth(4)
        h.Fit("fsum","LMER")

        if (fsum.GetNDF()>0):
            self.redChi2=fsum.GetChisquare()/fsum.GetNDF()
        else:
            self.redChi2=1000 
            #fsum.Draw("samel")


        return fsum



    
    def fitExp_cutoff (self, h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print "fit con gaus + exp+ cut-off"
        max = h.GetBinCenter( h.GetMaximumBin())
        convp=self.distConv
        print "distConv= ",self.distConv


        G0 = ROOT.TF1 ("G0","gaus",max-0.1,max+0.1)
        h.Fit("G0","LMER")
        G0.SetParameter(1,max)
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        G0.Draw("samel")

        
        
        #OK !!
        G1=ROOT.TF1 ("G1","gaus",convp-0.1,convp+0.1)
        G1.SetLineColor(4)
        h.Fit("G1","LMER")
        G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)

        f1= ROOT.TF1 ("f1","  (exp([0]+[1]*x))*(1./  (1+exp( (x-[2])/[3])  ))", convp+0.1,max+0.2)
        f1.SetParameters(2.3,4,max,0.01)

        f1.SetParLimits(2, max-sigma0/2., max+sigma0)
        #f1.FixParameter(2, mean0)
        h.Fit("f1","LMER")
        f1.SetLineColor(3)
        f1.Draw("samel")
        
        fsum=ROOT.TF1("fsum","gaus(0)+ (  ( exp([3]+[4]*x)*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  ))))",convp-0.2,max+0.2)

        #ok
        fsum.SetParameter(0, G1.GetParameter(0))
        #fsum.SetParameter(1, convp)
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))

        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)

        #ok
        #fsum.SetParLimits(0, 30,500)
        fsum.SetParLimits(0, 15,200)
       
        fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0) ) 
        fsum.SetParLimits(2, 0.02, sigma0) #!!!!!!!!!!!!!
        fsum.SetParLimits(5,  mean0-sigma0, mean0+sigma0)
        
        fsum.SetLineColor(2)
        fsum.SetLineWidth(4)
        h.Fit("fsum","LMER")
        if (fsum.GetNDF()>0):
            self.redChi2=fsum.GetChisquare()/fsum.GetNDF()
        else:
            self.redChi2=1000 


        return fsum

    def cerca_piccoAugerElectron (self, h, fLin, minX,par):
       # h e' l'istogramma in funzione della lunghezza lungo la traccia
       # restituisce le coord del picco auger 
       
    
        n_foundPeaks=0

        x_picco=0
        
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


        if (self.peakFinding==4): 
            # uso fit gauss + exp + cutoff
            self.fFit_histo=self.fitExp_cutoff(h)
            self.x_picco=self.fFit_histo.GetParameter(1)


            
        if (self.peakFinding==3):
            self.fFit_histo=self.fit2Gaussiane(h)
            self.x_picco=self.fFit_histo.GetParameter(1)



            
         # converto allo spazio xy    
                    
        x_lin= self.get_xDist(fLin, minX,self.x_picco)
        y_lin=par[0]*(x_lin**3)+par[1]*(x_lin**2)+(par[2]*x_lin)+par[3]
            
          
        newPoint=[x_lin,y_lin]
        return newPoint
         
      
  

