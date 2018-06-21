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
#import scipy
from scipy.interpolate import UnivariateSpline

import smoothing_passabassoSimo  as smooth_simo
from xpeSimo_ttree import *


import matplotlib.pyplot as plt
from gpdswswig.Recon import *
from gpdswswig.Geometry import *

from scipy.integrate import quad

#from numba import jit


class xpeSimo(object):

    def __init__(self,  raggioCut,dividiBins, baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw  ):

        
        self.track=-1 #!!!!!!!!!!!
        self.event_id=-100
        self.baricenter_X=-100
      
        self.baricenter_Y=-100
        self.conversion_point_X=-100
        self.conversion_point_Y=-100 
        self.phi0=-100
        self.phi1=-1090
        self.phiTang=-100
        
        
        self.xnew=0
        self.ynew=0 #  nuovo punto di conversione!!
               
              
        self.distConv=0.  #distanza  punto di conversione rec standard, lungo la traccia
        self.x_picco=0.   #
        self.distConvMC=0.#distanza  punto di conversione rec standard, lungo la traccia VERA (MC) 
        
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
        self.h1LCumul= ROOT.TH1F("h1LCumul","",int (self.nbinsX*10),self.x1,self.x2)

        
        self.h1L_smoothSimo= ROOT.TH1F("h1L_smoothSimo","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)
        self.h1L_ave= ROOT.TH1F("h1L_ave","",int (self.nbinsX/self.dividi_bins),self.x1,self.x2)

        self.h2 = ROOT.TH2F("h2","",self.nbinsX,self.x1,self.x2,self.nbinsY,self.y1,self.y2)
        self.profx=ROOT.TProfile("profx","profile",100,-2,2,"")


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
        self.f_spline3=ROOT.TF1()
        self.fLin2=ROOT.TF1("fLin2",self.f_dist, self.x1,self.x2,0)
        # self.f_splineScipy=ROOT.TF1()
        self.f_splineScipy=ROOT.TF1("f_splineScipy", self.eval_scipySpline_TF1, self.x1, self.x2,0 )
         
        #self.h2 = ROOT.TH2F()
        #self.gr1=ROOT.TGraphErrors() # aggiunto...
        #self.profx=ROOT.TProfile()
        self.newPoint=[0,0]
        self.gIon=ROOT.TGraph()
        self.MCconvPoint=ROOT.TMarker()

        """
        ROOT.SetOwnership(self.h2, False)
        ROOT.SetOwnership(self.gr1, False)
        ROOT.SetOwnership(self.profx, False)
        ROOT.SetOwnership(self.gIon, False)
        ROOT.SetOwnership(self.f_p3, False)
        ROOT.SetOwnership(self.f_spline3, False)
        ROOT.SetOwnership(self.fLin2, False)
        ROOT.SetOwnership(self.h1, False)
        ROOT.SetOwnership(self.h1L, False)
        ROOT.SetOwnership(self.h1L_ave, False)
        ROOT.SetOwnership(self.h1L_smoothSimo, False)
        ROOT.SetOwnership(self.outRootFile, False)
        ROOT.SetOwnership(self.f_splineScipy, False)
        """
        
        
        
        # parametri per rototraslazione!!!
        
        
        self.McInfo=-1

        # scipy spline:
        #self.uniSpline=UnivariateSpline
        #self.uniSplineDerivative= 0

        

    def fitScipy_spline(self,x,y,err):
        self.uniSpline= UnivariateSpline(x, y, w=err, s=10)
        self.uniSpilineDerivative=self.uniSpline.derivative()
        
        

            
    def eval_scipySpline_TF1(self, x):
        xx=float(x[0]) 
        return  self.uniSpline(xx) 

  
        
   
    
    def spline3(self,x,par):

        xx = float(x[0])
        xn = numpy.array( [par[0], par[1], par[2], par[3]], dtype=float)
        yn = numpy.array ([par[4], par[5], par[6], par[7]], dtype=float) 
        b1 =float(par[8])
        e1 =float( par[9])
        sp3=ROOT.TSpline3("sp3", xn, yn, 4, "b0e0", b1, e1);
        return sp3.Eval(xx)
        
    
    def f_dist(self,x):
        xx =float(x[0])
        #yy=math.sqrt( 1+( math.pow( self.f_p3.Derivative(xx),2)  ) )  #!!!!!!!!! fp3!!!!!    
        #yy=math.sqrt( 1+( math.pow( self.f_spline3.Derivative(xx),2)  ) )  #!!!!!!!!! fp3!!!!!    
        #yy=math.sqrt( 1+( math.pow( self.f_splineScipy.Derivative(xx),2)  ) )  #!!!!!!!!! spline con scipy!!!!!!!
        yy=math.sqrt( 1+( math.pow( self.uniSpilineDerivative(xx),2)  ) )  #!!!!!!!!! spline con scipy!!!!!!!
       
        return yy


    def f_distNew(self,x):
         yy=math.sqrt( 1+( math.pow( self.uniSpilineDerivative(x),2)  ) )  #!!!!!!!!! spline con scipy!!!!!!!
         return yy

    def distNew(self,x1,x2):

        dist=quad(self.f_distNew,x1,x2)
        return dist[0]
    

    
       
    def rec_simo(self):
        
        self.baricenter_X=self.track.barycenter().x()
        self.baricenter_Y=self.track.barycenter().y()
        self.conversion_point_X=self.track.absorptionPoint().x()
        self.conversion_point_Y=self.track.absorptionPoint().y()
        self.phi0=self.track.firstPassMomentsAnalysis().phi()
        self.phi1=self.track.secondPassMomentsAnalysis().phi()
        
        #if self.McInfo!=-1:
                    
        #creo liste con x,y,z di ogni pixel!!! forse meglio spostarlo fuori da questa classe e passare le liste!
        n_hits=self.track.numHits()
        hit=self.track.hits()
        #print ("n HITS = ",n_hits)
        
        x=numpy.array([0.]*n_hits)
        y=numpy.array([0.]*n_hits)
        adc=numpy.array([0.]*n_hits)
        
        for i in range (0,n_hits):
            x[i]=hit[i].x
            y[i]=hit[i].y
            adc[i]=hit[i].pulseHeight
            #print (x[i]," ",y[i]," adc = ",adc[i])

        
        
        xp,yp= self.rotoTraslate(x,y)

        
        phi0_bary=self.phi0-self.phi0
        phi1_bary=self.phi1 -self.phi0
               
                
        self.maxX= max(xp)
        self.minX=min(xp)
        self.maxY= max(yp)
        self.minY=min(yp)
                
       # self.h2 = ROOT.TH2F("h2","",self.nbinsX,self.x1,self.x2,self.nbinsY,self.y1,self.y2)
        #self.profx=ROOT.TProfile("profx","profile",100,-2,2,"")

        # ora e' importante!!!
        self.h1.Reset()
        self.h1L.Reset()
        self.h1L_smoothSimo.Reset()
        self.h1L_ave.Reset()
        self.h2.Reset()
        self.profx.Reset()

        
        
        #err=0.1*numpy.sqrt(adc)/adc # come normalizzo l'errore?????
        err= (adc/float(adc.max()))**0.5
        #print ("Err = ",err)
 
        
        for i in range (0,len(x)):
             self.h2.Fill(xp[i], yp[i],adc[i])
             self.profx.Fill(xp[i], yp[i],  numpy.sqrt(adc[i])) # peso con la radice dei conteggi... assumo che siano sempre >0!!! 
             self.h1.Fill(xp[i]-self.minX,adc[i])
        
        self.h2.GetXaxis().SetRangeUser(-0.8,1)
        self.h2.GetYaxis().SetRangeUser(-0.6, 0.6)        
        #self.gr1 = ROOT.TGraphErrors (len(xp),xp,yp,err,err )            
        #MC ionization track: !!!!! solo con files MC   FIXMEsimo!!!!
        #if self.McInfo!=-1:    
        xConvMC,yConvMC=self.rotoTraslate(self.McInfo.absorbtionPointX,self.McInfo.absorbtionPointY )
    
        self.MCconvPoint=ROOT.TMarker( xConvMC, yConvMC ,20)
        self.MCconvPoint.SetMarkerColor(6)
            
            
        ionXnp=numpy.array(self.McInfo.ionizationPosX)
        ionYnp=numpy.array(self.McInfo.ionizationPosY)
                   
        
        ion_xp,ion_yp=self.rotoTraslate(ionXnp, ionYnp)

        ionX=array('d',ion_xp)
        ionY=array('d',ion_yp)
        nIon=len(ionX)
        self.gIon=ROOT.TGraph(nIon,ionX,ionY) 

        
        
        
        #------------------ funzione solo con files MC!!!   
        
            
        #calcolo coodrdinate bary1 e 2 nel sistema roto traslato:       
        x_bary1,y_bary1=self.rotoTraslate(self.baricenter_X,self.baricenter_Y)
        self.bary1=ROOT.TMarker(x_bary1,y_bary1,20)
        self.bary1.SetMarkerColor(2)

        x_bary2,y_bary2=self.rotoTraslate(self.conversion_point_X,self.conversion_point_Y)
        self.bary2=ROOT.TMarker(x_bary2, y_bary2,20)
        self.bary2.SetMarkerColor(4)         
                
        m=(math.tan(phi0_bary))
        q=y_bary1-m*x_bary1

        m2=(math.tan(phi1_bary))
        q2=y_bary2-m2*x_bary2
       
        # fitto con spline scipy:
        err_fitMC=numpy.array([1.]*len(ion_xp))
        self.fitScipy_spline(ion_xp,ion_yp,err_fitMC)
        
        #compute X2:
        chi2=0.
        for i in range (0,len(ion_xp)):
            chi2=chi2+( (ion_yp[i]-self.uniSpline(ion_xp[i]) )**2)
        print ("======> chi2 = ",chi2)     
        
        
        if chi2>0.13:
            return -1
        
        #self.profx.Fit("f_p3","MERw")
                
        x_min_dist=0

        self.distConv=self.distNew(self.minX,0.) # in questo sitema di rif. il punto di conv e' in (0,0)
        self.distConvMC=self.distNew(self.minX, xConvMC) # in questo sitema di rif. il punto di conv e' in (0,0)
                
        for i in range (0,len(x)):
                    
             x_min_dist2=self.min_dist2(xp[i], yp[i], self.uniSpline, self.minX, self.maxX)
             y_min_dist2=self.uniSpline( x_min_dist2)
             distL2=self.distNew(self.minX,x_min_dist2 )
             radius=math.sqrt( (x_min_dist2-xp[i])**2 + (y_min_dist2-yp[i])**2)
             
             if radius <self.raggioCut:   # e se peso inversamente al raggio???
                # self.h1L.Fill(distL2,adc[i]/(2.7**(-radius**2)))
                # self.h1LCumul.Fill(distL2-self.distConvMC,adc[i]/(2.7**(-radius**2)))
                 self.h1L.Fill(distL2,adc[i])
                 self.h1LCumul.Fill(distL2-self.distConvMC,adc[i])

                 
 
        self.h1L_smoothSimo=self.smooth_simo(self.h1L)
        self.h1L_ave=self.media_mobile(self.h1L)
        
        
        # restituisce punto sulla traccia
        #self.newPoint= self.cerca_piccoAugerElectron(self.h1L_smoothSimo, fLin, minX,par)
        self.newPoint= self.cerca_piccoAugerElectron(self.h1L_smoothSimo, self.minX)
        #self.newPoint= self.cerca_piccoAugerElectron(self.h1L_ave, self.minX)
        #self.newPoint= self.cerca_piccoAugerElectron(self.h1L, self.minX)


        #self.newPoint= self.cerca_piccoAugerElectron(self.h1L, fLin2, minX)

        #m=3.*a*self.newPoint[0]**3+2.*b*self.newPoint[0]+c
        #m=self.f_p3.Derivative(self.newPoint[0])
        #print ("m=",m," m2 =",m2," diff = ",m-m2 )
                
        #self.phiTang=numpy.arctan(m)-phi       
                
        
        new_coord=self.undo_rotoTraslation (self.newPoint[0],self.newPoint[1])
     
        self.xnew=new_coord[0]
        self.ynew=new_coord[1]
        
                        

    #    del self.McInfo,  ionX,ion_xp,ionY,ion_yp
        
    def draw_simo(self):
        #============================================
        #    Draw!!!!!!!!!!
        #===========================================
        self.c_init.Clear()
        self.c_init.Divide(2,1)
        self.c_init.cd(1)
        
        self.h2.Draw("colZ")
        self.profx.Draw("samep")



       
        phi0_bary=self.phi0-self.phi0
        phi1_bary=self.phi1 -self.phi0
        #calcolo coodrdinate bary1 e 2 nel sistema roto traslato:       
        x_bary1,y_bary1=self.rotoTraslate(self.baricenter_X,self.baricenter_Y)
        self.bary1=ROOT.TMarker(x_bary1,y_bary1,20)
        self.bary1.SetMarkerColor(2)

        x_bary2,y_bary2=self.rotoTraslate(self.conversion_point_X,self.conversion_point_Y)
        self.bary2=ROOT.TMarker(x_bary2, y_bary2,20)
        self.bary2.SetMarkerColor(4)



        m=(math.tan(phi0_bary))
        q=y_bary1-m*x_bary1

        m2=(math.tan(phi1_bary))
        q2=y_bary2-m2*x_bary2
      
        #self.gr1.Draw("samep")

        #ellips=ROOT.TEllipse(self.baricenter_X,self.baricenter_Y, math.sqrt(self.mom2_long), math.sqrt(self.mom2_trans), 0,360, self.phi0*ROOT.TMath.RadToDeg())
        #ellips.SetFillStyle(0)
        #ellips.Draw()

        self.line1=ROOT.TF1("f1","[0]*x+[1]",self.minX, self.maxX)
        self.line1.SetParameters(m,q)

        self.line2=ROOT.TF1("f2","[0]*x+[1]",self.minX,self.maxX)
        self.line2.SetParameters(m2,q2)
        self.line2.SetLineColor(4)

        self.bary1.Draw()   
        self.bary2.Draw()
           
        self.line1.Draw("samel")
        self.line2.Draw("samel")


        xFit=numpy.linspace(self.minX, self.maxX,1000)
        yFit=self.uniSpline(xFit)
        
        gfit=ROOT.TGraph(1000,array('f',xFit),array('f',yFit))
        gfit.SetMarkerColor(6)
        gfit.SetLineColor(6)
        gfit.SetLineWidth(2)
        gfit.Draw("pl")
        
        #self.f_splineScipy.Draw("samel") !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        PeakImpP=ROOT.TMarker( self.newPoint[0],  self.newPoint[1],34)
        PeakImpP.SetMarkerColor(6)
        PeakImpP.SetMarkerSize(2.1)
        PeakImpP.Draw()   
                              
        self.h1L.GetXaxis().SetRangeUser(-0.4,1.5)
        self.h1L.SetLineColor(1)
        self.h1L.SetLineWidth(1)
            
        # conv point MC:
        if self.McInfo!=-1:
            self.MCconvPoint.Draw()
            self.gIon.Draw("*")
                       
        

        
        ###########
        # pad2
        self.c_init.cd(2)
        self.h1L.Draw("hist")
        convPoint = ROOT.TMarker(self.distConv,0,22)
        convPoint.SetMarkerColor(4)
        convPoint.Draw()

        convPointMC = ROOT.TMarker(self.distConvMC,0,22)
        #convPoint.SetMarkerColor(2)
        convPointMC.Draw()

        peakPoint = ROOT.TMarker(self.x_picco,0,22)
        peakPoint.SetMarkerColor(2)
        peakPoint.Draw()
            
        
        
        self.h1L_smoothSimo.SetLineColor(4)
        self.h1L_smoothSimo.SetLineWidth(2)
        self.h1L_smoothSimo.Draw("sames")
            
            
        self.h1L_ave.SetLineColor(8)
        self.h1L_ave.SetLineWidth(4)
        #self.h1L_ave.Draw("sames")

        if self.peakFinding>1 and self.peakFinding!=5:
            self.fFit_histo.Draw("samel")

        #if self.peakFinding==2:
        #    self.fFit_histo.Draw("samel")   

        #if self.peakFinding==2:
        #    self.fFit_histo.Draw("samel")   
             
        self.c_init.Update() #!!!!!!!!!!!!!!!!!!!!!!!!!
        self.outRootFile.cd()
        self.c_init.Write()
            
        valore = input('continue?')
                
        
        return 1


          
    #@jit     
    def min_dist2 (self, pixel_x, pixel_y, func, xmin,xmax): 
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
                #y=func.Eval(x)
                y=func(x)
                dist=math.sqrt( (x-pixel_x)*(x-pixel_x)+(y-pixel_y)*(y-pixel_y))
                if (dist<=min): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aggiunto <= !!!
                    min=dist
                    x_min=x

            x_start=x_min-dx1
            x_stop=x_min+dx1
        #end for n_iters

        #print ("min dist = ",min, "xmin = ",x_min)
        return x_min


   
    #jit    
    def get_xDist (self,x_min, l):
        """ cerco l'estremo superiore tale che l'integrale tra xmin e b sia ==l   
             SAREBBE DA CALCOLARE ANALITICAMENTE!!!!!!!, ma per ora mi accontento di calcolo numerico... 
        """
        x_max=l
        n_steps=30.
        dx=(x_max-x_min)/n_steps

        sup_min=x_min
        xfound=0
        
        for n_iter in range (0,2):

            min_d=1000
            for i in range (0, int( n_steps)):
                sup=sup_min+i*dx 
                d_calc=  math.fabs(self.distNew(x_min,sup)-l)
                
                
                if (d_calc<min_d):
                    min_d=d_calc
                    xfound=sup
            
            # setta i parametri per la prossima iterazione...           
            sup_min=xfound-dx
            sup_max=xfound+dx
            dx=(sup_max-sup_min)/n_steps        
            
            
        #print ("returning xfound = ",xfound)
        return xfound 

    

   # @jit
    def smooth_simo (self, h1):
       # ficco i valore dellh su un vettore
       nbins=h1.GetNbinsX()
       x=numpy.array([0.]*nbins, dtype=float)
       for i in range (1,nbins):
           x[i]=h1.GetBinContent(i)

       #print ("creo oggetto...")
       simoFilter=smooth_simo.filtra_segnale()

       #print ("init filtro.. ")
       aa=simoFilter.init_filtro();
       y=simoFilter.filtra(x)
       hF=h1.Clone()
       hF.Reset()
       for i in range (1,nbins):
           hF.SetBinContent(i,y[i])
           #print ("filtered: i= ",i," y = ",y[i])

       return hF

    def media_mobile (self, h1):        
       nbins=h1.GetNbinsX()
       hF=h1.Clone()
       hF.Reset()
       for i in range (2,nbins-2):
             
           #ave=  (h1.GetBinContent(i)+ h1.GetBinContent(i-1)+h1.GetBinContent(i+1))/3.
           ave=  (h1.GetBinContent(i)+ h1.GetBinContent(i-1)+h1.GetBinContent(i+1)+h1.GetBinContent(i-2)+h1.GetBinContent(i+2)  )/5.
           hF.SetBinContent(i,ave)
       return hF      

    #@jit
    def undo_rotoTraslation (self,xp,yp):

       phi=self.phi0
       x0=self.conversion_point_X
       y0=self.conversion_point_Y
       
       dx=xp*numpy.sin(phi)+yp*numpy.cos(phi)
       dy=xp*numpy.cos(phi)-yp*numpy.sin(phi)

       x=dx+x0  
       y=dy+y0
       new_coord=[x,y]
       return new_coord        
   


    def firstBinAboveTh(self,h1):
       th=50.
       x_max=0 
       for i in range(1,h1.GetNbinsX()):
           #print ("i=",i, "cont = ",h1.GetBinContent(i))
           if h1.GetBinContent(i)>th:
               x_max=h1.GetBinCenter(i)
               #print (" trovato i=",i,"x_max=",x_max)
               break

       G0 = ROOT.TF1 ("G0","gaus",0.2,0.6)   
       return x_max
   
    def fit2Gaussiane (self, h1):

       G0 = ROOT.TF1 ("G0","gaus",0.2,0.6)
       h1.Fit("G0","MER")
       mean0=G0.GetParameter(1)
       sigma0=G0.GetParameter(2)

       G2 = ROOT.TF1 ("G2","gaus",0.25,0.5)
       G2.SetParameter(0, 9.74680e+02 )
       G2.SetParameter(1,  3.17394e-01)
       G2.SetParameter(2, 5.47369e-02 )

      # G2.SetParLimits(1, mean0-sigma0/2., mean0+sigma0/2.)
       #G2.FixParameter(1)
       #G2.FixParameter(2)
       
       h1.Fit("G2","MER")
        

       G1 = ROOT.TF1 ("G1","gaus",0.05,0.35)
       G1.SetParameter(0,3.33294e+02 )
       G1.SetParameter(1, 1.35906e-01 )
       G1.SetParameter(2, 5.30043e-02  )

      # G1.SetParLimits(1, 0.03,0.25)
      # G1.SetParLimits(2, 0.02,0.1)
       
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

       #Gsum.SetParLimits(1, 0.03,0.3)
       #Gsum.SetParLimits(2, 0.02,0.1)
       #Gsum.SetParLimits(4,  mean0-sigma0/2., mean0+sigma0/2.)
       
            
       Gsum.SetLineColor(2)
       h1.Fit("Gsum","MER")
       
        
       return Gsum
       


    def fitGaus_cutoff (self, h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        #print ("fit con cut-off")
        
        max = h.GetBinCenter( h.GetMaximumBin())
        convp=self.distConv
        #print ("distConv= ",self.distConv)

        
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
        
        print ("fit con gaus + exp+ cut-off")
        max = h.GetBinCenter( h.GetMaximumBin())
        convp=self.distConv
        #print ("distConv= ",self.distConv)


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



        #print (" sigma G1 = ",fsum.GetParameter(2)  )
            
        return fsum


    
    def cerca_piccoAugerElectron (self, h, minX):
       # h e' l'istogramma in funzione della lunghezza lungo la traccia
      # restituisce le coord del picco auger 
       
    
        n_foundPeaks=0

        x_picco=0
        
        #cerca picchi!
        if (self.peakFinding==1):
            #print (" inizio ROOT.Tspectrum")
            s=ROOT.TSpectrum(self.maxnP)
            n_foundPeaks=s.Search(h,self.Psigma,"",self.Pthr)               
            x_peaks_raw=s.GetPositionX()
            x_peaks=[]
        
            for i in range (0,n_foundPeaks):
                x_peaks.append(x_peaks_raw[i])

            x_peaks.sort()

            #for i in range (0,n_foundPeaks):
            #   print ("x peak found = ",x_peaks[i] )
        
            if n_foundPeaks >=2:
                self.x_picco=x_peaks[0]
            else:
                newPoint=[-100,100]
                return newPoint
        #------------------------        
                
        if (self.peakFinding==2): 
            # uso fit doppia gaussiana
            #self.fFit_histo=self.fitGaus_cutoff(h)
            self.fFit_histo=self. fit2Gaussiane(h)
           
            self.x_picco=self.fFit_histo.GetParameter(1)


        if (self.peakFinding==4): 
            # uso fit gauss + exp + cutoff
            self.fFit_histo=self.fitExp_cutoff(h)
            self.x_picco=self.fFit_histo.GetParameter(1)


            
        if (self.peakFinding==3):
            self.fFit_histo=self.fit2Gaussiane(h)
            self.x_picco=self.fFit_histo.GetParameter(1)


        if (self.peakFinding==5):
            #self.fFit_histo=self.fit2Gaussiane(h)
            self.x_picco=self.firstBinAboveTh(h)

            
         # converto allo spazio xy    
                    
        x_lin= self.get_xDist(minX,self.x_picco)
        #y_lin=par[0]*(x_lin**3)+par[1]*(x_lin**2)+(par[2]*x_lin)+par[3]
        #y_lin=self.f_spline3.Eval(x_lin)
        y_lin=self.f_splineScipy.Eval(x_lin)
        
        #print ("x_picco = ",self.x_picco, "x_lin = ",x_lin)
        
        #print ("===============>>>>>>>>>>>> y_lin=",y_lin,"  y_lin2=",y_lin2,"  diff= ",y_lin-y_lin2)
          
        newPoint=[x_lin,y_lin]
        return newPoint
         
      
  
    
    def rotoTraslate(self, x,y):
          
        x0=self.conversion_point_X
        y0=self.conversion_point_Y
        phi=self.phi0
           
        dx = (x - x0)
        dy = (y - y0)
        xp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
        yp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy
        #aaa=[xp,yp]
        return xp,yp
