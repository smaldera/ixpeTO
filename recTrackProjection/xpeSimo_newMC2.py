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
#from fit_distro_simo  import fitDistrib  as fitSimo
import fit_distro_simo  

import matplotlib.pyplot as plt
from gpdswswig.Recon import *
from gpdswswig.Geometry import *
from gpdswswig.MonteCarlo import *


from numba import jit


class xpeSimo(object):
    def __del__(self):
        print ("destructor")
        self.deleteHistos()


    def __init__(self, track,   raggioCut,dividiBins, baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw  ):
   # def __init__(self,  raggioCut,dividiBins, baryPadding, findMaxAlg, pcubo, maxnP, Psigma, Pthr, draw  ):

       
        
        self.event_id=-1
        self.baricenter_X=track.barycenter().x()
        print ("bary X = ",self.baricenter_X)
        self.baricenter_Y=track.barycenter().y()
        self.conversion_point_X=track.absorptionPoint().x()
        self.conversion_point_Y=track.absorptionPoint().y()
       
        self.phi0=track.firstPassMomentsAnalysis().phi()
        self.phi1=track.secondPassMomentsAnalysis().phi()
        self.phiTang=0
        self.track=track #!!!!!!!!!!!
    
  
        
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
        self.x2=4 
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
        self.h2 = ROOT.TH2F("h2","",self.nbinsX,self.x1,self.x2,self.nbinsY,self.y1,self.y2)    
        print ("INIT, nbinsX=",self.nbinsX)
        

        
        # funzione fit istogramma long. cariche
        self.fFit_histo=0
        self.redChi2=0
        self.sumpars=array('d',[0.]*7) 
        
        #self.c_init=0
        self.c_init= ROOT.TCanvas("cc","cc", 2000,1000) 
        self.c_init.Divide(2,1)
        self.c_init.Draw()
        self.outRootFile=ROOT.TFile()

        
        self.f_p3=ROOT.TF1()
        self.f_spline3=ROOT.TF1()
        self.fLin2=ROOT.TF1()
        self.minX=0
        self.maxX=0
        
        self.fFit=ROOT.TF1() # funzione di fit della traccia da usare dappertutto!!!!!
        
        self.gr1=ROOT.TGraphErrors() # aggiunto...
        self.newPoint=[0,0]

        # parametri per rototraslazione!!!
        
        self.phi=self.phi0
        self.x0=self.conversion_point_X
        self.y0=self.conversion_point_Y

        self.McInfo=-1


        # scipy spline:
        self.uniSpiline=0

    def fitScipy_spline(self,x,y,err):
            self.uniSpline= UnivariateSpline(x, y, w=err, k=5,   s=None)
            
            
    def eval_scipySpline_TF1(self, x):
            xx=float(x[0]) 
            return  self.uniSpline(xx) 

  
        
    """
    def fit_spline(self, zero_suppression=12):
        #To be moved into recon.
        
        from scipy.interpolate import UnivariateSpline
        _mask = self.cluster.adc_values >= zero_suppression
        x = self.cluster.x[_mask]
        y = self.cluster.y[_mask]
        adc_values = self.cluster.adc_values[_mask]
        weights = (adc_values/float(adc_values.max()))**0.5
        dx = (x - self.baricenter_X)
        dy = (y - self.baricenter_Y)
        xp = numpy.cos(self.phi0)*dx + numpy.sin(self.phi0)*dy
        yp = -numpy.sin(self.phi0)*dx + numpy.cos(self.phi0)*dy
        #s = UnivariateSpline(xp, yp, w=weights, s=0.5)
        s = UnivariateSpline(xp, yp, w=weights, s=1)
       
        _xp = numpy.linspace(xp.min(), xp.max(), 25)
        _yp = s(_xp)
        dx = numpy.cos(-self.phi0)*_xp + numpy.sin(-self.phi0)*_yp
        dy = -numpy.sin(-self.phi0)*_xp + numpy.cos(-self.phi0)*_yp
        x = dx + self.baricenter_X
        y = dy + self.baricenter_Y 
        plt.plot(x, y, '-', lw=2, color='black')
    """

    
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
        yy=math.sqrt( 1+( math.pow( self.fFit.Derivative(xx),2)  ) )      
        return yy

    def deleteHistos(self):
        del  self.h1      
        del  self.h1L
        del  self.h1L_smoothSimo
        del  self.h1L_ave
        del  self.h2
        
    def ResetHistos(self):     
        self.h1.Reset()
        self.h1L.Reset()
        self.h1L_smoothSimo.Reset()
        self.h1L_ave.Reset()
        #self.h2.Reset()


    def creaFunzioniFit(self):     

        self.f_p3=ROOT.TF1("f_p3","[0]*x*x*x+[1]*x*x+[2]*x+ ([3]- [0]*[4]*[4]*[4]-[1]*[4]*[4]-[2]*[4] )", self.minX, self.maxX) # pol3  per bary2
       # self.f_p3=ROOT.TF1("f_p3"," ( [0]*x*x*x+[1]*x*x+[2]*x+ ([3]- [0]*[4]*[4]*[4]-[1]*[4]*[4]-[2]*[4] ))*(x>0)+0*(x<=0)", minX, maxX) # pol3  per bary2
        #self.f_p3=ROOT.TF1("f_p3","[0]*x+[1]", minX, maxX) # pol3  per bary2
               
        # fissa il passaggio per il punto di conv. della rec standard
        #self.f_p3.FixParameter(3,y_bary2)
        #self.f_p3.FixParameter(4,x_bary2)

        # fissa il passaggio per il punto di conv. MC
        self.f_p3.FixParameter(3,self.yConvMC)
        self.f_p3.FixParameter(4,self.xConvMC)

        
        #self.f_p3.SetParLimits(3,y_bary2-self.baryPadding ,y_bary2+self.baryPadding)
        #self.f_p3.SetParLimits(4,x_bary2-self.baryPadding,x_bary2+self.baryPadding)

        #if self.pcubo==0:
        #    self.f_p3.FixParameter(0,0) # il cubo diventa una parabola!
        #    self.f_p3.FixParameter(1,0) # il cubo diventa una retta!!!!!
        
        ########### SPLINE 4 NODI ROOT
        #uso la spline3(4nodi)!!!!!!!!!!!!!!!!!!!!!!
        self.f_spline3=ROOT.TF1("f_spline3",self.spline3, self.minX, self.maxX, 10) # spline3  per ora nessun controllo su passagio per/vicino  p.conv.

        deltaX=(self.maxX-self.minX)/5.
        

        self.f_spline3.SetParameters(1.*deltaX ,2.*deltaX,3.*deltaX,4.*deltaX,self.f_p3.Eval(1.*deltaX ),self.f_p3.Eval(2.*deltaX),self.f_p3.Eval(3.*deltaX), self.f_p3.Eval(4.*deltaX),self.f_p3.Derivative(self.minX),self.f_p3.Derivative(self.maxX))#parametri iniziali!!!!!
        

        
        #self.f_spline3=ROOT.TF1("f_spline3","[0]*x*x*x+[1]*x*x+[2]*x+ ([3]- [0]*[4]*[4]*[4]-[1]*[4]*[4]-[2]*[4] )", self.minX, self.maxX) # pol3  per bary2 ?????????
        #self.f_spline3.SetParLimits(3,y_bary2-self.baryPadding ,y_bary2+self.baryPadding)
        #self.f_spline3.SetParLimits(4,x_bary2-self.baryPadding,x_bary2+self.baryPadding)
        self.f_spline3.SetLineColor(6);
        self.f_spline3.SetLineWidth(5);

        ###############
        # SCIPY UNIVARIATE SPLINE
        
    """       
    def m_p3analitico(self):
        par3=self.f_p3.GetParameter(3)
        par4=self.f_p3.GetParameter(4)
        a=self.f_p3.GetParameter(0)
        b=self.f_p3.GetParameter(1)
        c= self.f_p3.GetParameter(2)
        d= par3-a*par4**3-b*par3**2-c*par4    
        m=3.*a*self.x_bary2**3+2.*b*self.x_bary2+c
        return m
    """ 


        
    def rec_simo(self):

        print ("")
        print ("==================== ")
        print ("")
        
        import time
        import math
       
        if self.McInfo!=-1:
            print ("Xmc=",self.McInfo.absorbtionPointX)
            
        
        print("MC_energy=",self.McInfo.photonEnergy)
        print("MC_PEenergy=",self.McInfo.photoElectronEnergy)
        print("MC_PEphi=",self.McInfo.photoElectronPhi)
        print("MC_PEtheta=",self.McInfo.photoElectronTheta)
        print("MC_AugerEnergy=",self.McInfo.augerEnergy)

        
        #creo liste con x,y,z di ogni pixel!!! forse meglio spostarlo fuori da questa classe e passare le liste!
        n_hits=self.track.numHits()
        hit=self.track.hits()
        
        x=numpy.array([0.]*n_hits)
        y=numpy.array([0.]*n_hits)
        adc1=numpy.array([0.]*n_hits)
        
        for i in range (0,n_hits):
            x[i]=hit[i].x
            y[i]=hit[i].y
            adc1[i]=hit[i].pulseHeight
                   
        
        
        #self.x0=self.McInfo.absorbtionPointX
        #self.y0= self.McInfo.absorbtionPointY
        #self.phi=self.McInfo.photoElectronPhi

        xAbs_rot,yAbs_rot=self.rotoTraslate(self.McInfo.absorbtionPointX  ,  self.McInfo.absorbtionPointY )
        if xAbs_rot>0:
            self.phi=self.phi-3.1415
        
        
        xp1,yp1= self.rotoTraslate(x,y)
        xp,yp,adc=self.ordinateByX(xp1,yp1,adc1)

                    
        phi0_bary=self.phi0-self.phi
        phi1_bary=self.phi1 -self.phi
                       
        self.maxX= max(xp)
        self.minX=min(xp)
        self.maxY= max(yp)
        self.minY=min(yp)
        self.ResetHistos()  # probabilmente inutile, giusto per essere sicuro che siano vuoti
       
                
        #err=0.1*numpy.sqrt(adc)/adc # come normalizzo l'errore?????
        err= (adc/float(adc.max()))**0.5       
        
        for i in range (0,len(x)):
             self.h2.Fill(xp[i], yp[i],adc[i])
             self.h1.Fill(xp[i]-self.minX,adc[i])
        
        
       
        ionXnp=numpy.array(self.McInfo.ionizationPosX)
        ionYnp=numpy.array(self.McInfo.ionizationPosY)
        ionZnp=numpy.array(self.McInfo.ionizationPosZ)
         
        ion_xp1,ion_yp1=self.rotoTraslate(ionXnp, ionYnp)
        ion_xp,ion_yp,ion_zp=self.ordinateByX(ion_xp1,ion_yp1,ionZnp)

        ionX=array('d',ion_xp)
        ionY=array('d',ion_yp)
        ionZ=array('d',ion_zp)
        nIon=len(ionX)

        x_all1= numpy.concatenate((ion_xp1,xp1),axis=None )
        y_all1=numpy.concatenate((ion_yp1,yp1),axis=None )
        w_all1=numpy.concatenate( (numpy.array([1]*len(ion_xp1)), ((adc1/float(adc1.max()))**0.5)),axis=None  )
        x_all, y_all,w_all=self.ordinateByX(x_all1,y_all1,w_all1)
     
        

        self.gr1 = ROOT.TGraphErrors (len(xp),xp,yp,err,err )
        self.gIon=ROOT.TGraph(nIon,ionX,ionY) 
        self.gIon_xz=ROOT.TGraph(nIon,ionX,ionZ)
        self.gIon_yz=ROOT.TGraph(nIon,ionY,ionZ)
        
        # fit angolo theta:
        rettaZ=ROOT.TF1("rettaZ","[0]+[1]*x",-1,1)
        self.gIon_xz.Fit("rettaZ","MER")
        thetaSimo=numpy.arctan(rettaZ.GetParameter(0))
        print("thetaSimo = ",ROOT.TMath.RadToDeg()*thetaSimo)
        
        minX_ion= min(ionX)
        maxX_ion= max(ionX)
        
        #calcolo coodrdinate bary1 e 2 nel sistema roto traslato:       
        self.x_bary1,self.y_bary1=self.rotoTraslate(self.baricenter_X,self.baricenter_Y)
        self.x_bary2, self.y_bary2=self.rotoTraslate(self.conversion_point_X,self.conversion_point_Y)
        self.xConvMC,self.yConvMC=self.rotoTraslate(self.McInfo.absorbtionPointX,self.McInfo.absorbtionPointY )  


       
            
        #self.gr1.Fit("f_spline3","MRw") ##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #self.gIon.Fit("f_spline3","MRw")
        #self.gIon.Fit("f_p3","MRw")       

        # fitto con spline scipy:      
        #ionZ=array('d',[0.1]*len(ionX) )
        #self.fitScipy_spline( ionX,ionY,ionZ)  
        self.fitScipy_spline( x_all,y_all,w_all)
        # creo TF1 della spline: (questo lo devo lasciare qua, esiste solo con il fit!!!!!)
        self.f_splineScipy=ROOT.TF1("f_splineScipy", self.eval_scipySpline_TF1, self.minX,self.maxX,0 ) # spline3  per ora nessun controllo su passagio per/vicino  p.conv.
                 
                
        
        self.creaFunzioniFit()      
        #self.fFit=self.f_splineScipy
        self.fFit=self.f_p3
        self.gIon.Fit("f_p3","ME") 
        
        fLin2=ROOT.TF1("fLin2",self.f_dist, self.minX, self.maxX,0)
        
        self.distBary=fLin2.Integral(self.minX, self.x_bary1)  
         
        
        #distCorrection=math.fabs( 1./math.cos(thetaSimo-3.1415/2.))
        distCorrection=1
        
       
        print("distCorrection = ",distCorrection, "theta = ",self.McInfo.photoElectronTheta*ROOT.TMath.RadToDeg(), " thetaSimo = ", thetaSimo*ROOT.TMath.RadToDeg() )
        r0=0.7
        for i in range (0,len(xp)):
                         
             
             x_min_dist2=self.min_dist2(xp[i], yp[i], self.fFit, self.minX, self.maxX)
             y_min_dist2=self.fFit.Eval( x_min_dist2)
            
             distL2=(fLin2.Integral(self.minX, x_min_dist2)-    self.distBary)   *distCorrection     # d=o a dist baricentro    
             radius=math.sqrt( (x_min_dist2-xp[i])**2 + (y_min_dist2-yp[i])**2)
             
             # if radius <self.raggioCut:   # e se peso inversamente al raggio???
             #     self.h1L.Fill(distL2,adc[i]) 
             
             self.h1L.Fill(distL2,adc[i]*math.exp(-radius/r0)) 
             # self.h1L.Fill(distL2,adc[i]) 



        ################### prova spot
        """
        minIon_d=fLin2.Integral(self.minX, xp[0])
        maxIon_d=fLin2.Integral(self.minX, max(xp))

        
        print("minIon_d=",minIon_d,"  max = ",maxIon_d)
        
        d_i=numpy.linspace(minIon_d,maxIon_d,200)
        print ("D_i=",d_i)

      
        
        for i in range (0,len(d_i)):
             x_i= self.get_xDist(fLin2,self.minX, d_i[i])  
             distIon=(d_i[i]-self.distBary)*distCorrection 
             y_i=self.fFit.Eval(x_i)                     
             for jj in range (0,len(xp)):
                radiusIon=math.sqrt( (x_i-xp[jj])**2 + (y_i-yp[jj])**2)
               # self.h1L.Fill(distIon,adc[jj]*math.exp(-radiusIon/r0))
        """    


            


        ####################





            
        self.distConv=fLin2.Integral(self.minX, self.x_bary2)-self.distBary 
        self.distConvMC=fLin2.Integral(self.minX, self.xConvMC)-self.distBary
        self.h1L_smoothSimo=self.smooth_simo(self.h1L)
        self.h1L_ave=self.media_mobile(self.h1L)
        self.newPoint= self.cerca_piccoAugerElectron(self.h1L_smoothSimo, fLin2, self.minX)
        

      
        m=self.fFit.Derivative(self.newPoint[0])  
        self.phiTang=numpy.arctan(m)-self.phi       
                
        
        new_coord=self.undo_rotoTraslation (self.newPoint[0],self.newPoint[1])
        self.xnew=new_coord[0]
        self.ynew=new_coord[1]
        
        
                        
        

        
    def draw_simo(self):
        #============================================
        #    Draw!!!!!!!!!!
        #===========================================
        self.c_init.Clear()
       # self.c_init.Divide(2,1)
        self.c_init.Divide(2,2)
       
        self.c_init.cd(1)

        self.h2.GetXaxis().SetRangeUser(-0.8,1)
        self.h2.GetYaxis().SetRangeUser(-0.6, 0.6)
        self.h2.Draw("colZ")
        #self.gr1.Draw("samep")

        #ellips=ROOT.TEllipse(self.baricenter_X,self.baricenter_Y, math.sqrt(self.mom2_long), math.sqrt(self.mom2_trans), 0,360, self.phi0*ROOT.TMath.RadToDeg())
        #ellips.SetFillStyle(0)
        #ellips.Draw()
        bary1Marker=ROOT.TMarker(self.x_bary1,self.y_bary1,20)
        bary2Marker=ROOT.TMarker(self.x_bary2,self.y_bary2,20)
      
        bary2Marker.SetMarkerColor(4)
        bary1Marker.SetMarkerColor(2)
        bary1Marker.Draw()   
        bary2Marker.Draw()



        phi0_bary=self.phi0-self.phi
        phi1_bary=self.phi1 -self.phi
        m=(math.tan(phi0_bary))
        q=self.y_bary1-m*self.x_bary1
        m2=(math.tan(phi1_bary))
        q2=self.y_bary2-m2*self.x_bary2
             
        line1=ROOT.TF1("f1","[0]*x+[1]",self.minX, self.maxX)
        line1.SetParameters(m,q)

        line2=ROOT.TF1("f2","[0]*x+[1]",self.minX,self.maxX)
        line2.SetParameters(m2,q2)
        line2.SetLineColor(4)


        
        line1.Draw("samel")
        line2.Draw("samel")       
        self.fFit.Draw("samel")
        
        

        PeakImpP=ROOT.TMarker( self.newPoint[0],  self.newPoint[1],34)
        PeakImpP.SetMarkerColor(6)
        PeakImpP.SetMarkerSize(2.1)
        PeakImpP.Draw()   
                              
        self.h1L.GetXaxis().SetRangeUser(-1,2)

        
        self.h1L.SetLineColor(1)
        self.h1L.SetLineWidth(1)
            
        # conv point MC:
        if self.McInfo!=-1:
            
 
            MCconvPoint=ROOT.TMarker( self.xConvMC, self.yConvMC ,48)
            MCconvPoint.SetMarkerColor(6)
            MCconvPoint.Draw()
            self.gIon.Draw("*")         


        
        ###########
        # pad2
        
        self.c_init.cd(2)
        self.h1L.Draw("hist")
        
        convPoint = ROOT.TMarker(self.distConv,0,22)
        convPoint.SetMarkerColor(2)
        convPoint.Draw()

        convPointMC = ROOT.TMarker(self.distConvMC,0,20)
        convPointMC.SetMarkerColor(6)
        convPointMC.Draw()

        
        peakPoint = ROOT.TMarker(self.x_picco-self.distBary,0,22)
        peakPoint.SetMarkerColor(4)
        peakPoint.Draw()
        
        
        self.h1L_smoothSimo.SetLineColor(4)
        self.h1L_smoothSimo.SetLineWidth(2)
        self.h1L_smoothSimo.Draw("sames")
            
            
        #self.h1L_ave.SetLineColor(8)
        #self.h1L_ave.SetLineWidth(4)
        #self.h1L_ave.Draw("sames")

        if self.peakFinding==3:
            self.fFit_histo.Draw("samel")

        if self.peakFinding==2:
            self.fFit_histo.Draw("samel")   

        ###########
        # pad3
        self.c_init.cd(3)
        self.gIon_xz.Draw("a*")
        ###########
        # pad4
        self.c_init.cd(4)
        self.gIon_yz.Draw("a*")

        






            
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
                y=func.Eval(x)
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
    def get_xDist (self,fLin2, x_min, l):
        """ cerco l'estremo superiore tale che l'integrale tra xmin e b sia ==l   
             SAREBBE DA CALCOLARE ANALITICAMENTE!!!!!!!, ma per ora mi accontento di calcolo numerico... 
        """
        x_max=l
        n_steps=30.
        dx=(x_max-x_min)/n_steps

        sup_min=x_min
        xfound=0

        print ("x_min=",x_min)
        
        for n_iter in range (0,2):

            min_d=1000
            for i in range (0, int( n_steps)):
                sup=sup_min+i*dx 
                d_calc=  math.fabs(fLin2.Integral(x_min,sup)-l)
            
                if (d_calc<min_d):
                    min_d=d_calc
                    xfound=sup
            
            # setta i parametri per la prossima iterazione...           
            sup_min=xfound-dx
            sup_max=xfound+dx
            dx=(sup_max-sup_min)/n_steps        
            
            
        print ("returning xfound = ",xfound)
        return xfound 

    

    @jit
    def smooth_simo (self, h1):
       # ficco i valore dellh su un vettore
       nbins=h1.GetNbinsX()
       x=numpy.array([0.]*nbins, dtype=float)
       for i in range (1,nbins):
           x[i]=h1.GetBinContent(i)

       print ("creo oggetto...")
       simoFilter=smooth_simo.filtra_segnale()

       print ("init filtro.. ")
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
             
             ave=  (h1.GetBinContent(i)+ h1.GetBinContent(i-1)+h1.GetBinContent(i+1))/3.
             hF.SetBinContent(i,ave)
       return hF      

    @jit
    def undo_rotoTraslation (self,xp,yp):

       phi=self.phi
       x0=self.x0
       y0=self.y0
       
       dx=xp*numpy.sin(phi)+yp*numpy.cos(phi)
       dy=xp*numpy.cos(phi)-yp*numpy.sin(phi)

       x=dx+x0  
       y=dy+y0
       new_coord=[x,y]
       return new_coord        

   
    
    def cerca_piccoAugerElectron (self, h, fLin2, minX):
       # h e' l'istogramma in funzione della lunghezza lungo la traccia
      # restituisce le coord del picco auger 
       
    
        n_foundPeaks=0

        x_picco=0
        
        #cerca picchi!
        if (self.peakFinding==1):
            print (" inizio ROOT.Tspectrum")
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
                self.x_picco=x_peaks[0]  +self.distBary      #!!!!!!!!!!!1
            else:
                newPoint=[-100,100]
                return newPoint
        #------------------------        
        fitSimo= fit_distro_simo.fitDistrib()       
        if (self.peakFinding==2):

            #fitSimo= fit_distro_simo.fitDistrib()
            
            # uso fit doppia gaussiana
            
            #self.fFit_histo=fitSimo.fit2Gaussiane(h)
            self.fFit_histo=fitSimo.super_Fit2Gaussiane(h)
            self.x_picco=self.fFit_histo.GetParameter(1)  +self.distBary      #!!!!!!!!!!!1


        if (self.peakFinding==4): 
            # uso fit gauss + exp + cutoff
            #self.fFit_histo=self.fitExp_cutoff(h)
            self.fFit_histo=fitSimo.fitExp_cutoff(h)
            self.x_picco=self.fFit_histo.GetParameter(2)  +self.distBary      #!!!!!!!!!!!1


            
        if (self.peakFinding==3):
            self.fFit_histo=fitSimo.fitExpCutoff_gaus(h)
            self.x_picco=self.fFit_histo.GetParameter(1)  +self.distBary      #!!!!!!!!!!!1


        if (self.peakFinding==5):
            print("cerca picco simo")
            self.x_picco=fitSimo.myPeakFind(h)  +self.distBary      #!!!!!!!!!!!1
            

            
         # converto allo spazio xy    
   
        x_lin= self.get_xDist(fLin2, minX,self.x_picco)
        y_lin=self.fFit.Eval(x_lin)
        
        print ("x_picco = ",self.x_picco, "x_lin = ",x_lin)
        
        #print ("===============>>>>>>>>>>>> y_lin=",y_lin,"  y_lin2=",y_lin2,"  diff= ",y_lin-y_lin2)
          
        newPoint=[x_lin,y_lin]
        return newPoint
         


    def ordinateByX(self,x,y,z):
    
       
       from operator import itemgetter
       dict_value={}
       for i in range (0,len(x)):
           dict_value[i]=x[i]
           mySort=sorted(dict_value.items(), key=itemgetter(1))

       index_sorted2= [val[0]  for  val in mySort]
       x_sorted2= numpy.array([val[1]  for  val in mySort])
       #y_sorted2=numpy.array([0.]*len(x),float)
       #z_sorted2=numpy.array([0.]*len(x),float)


       y_sorted2=numpy.array( [y[i]  for i in index_sorted2])
       z_sorted2=numpy.array( [z[i]  for i in index_sorted2])

       #print(index_sorted2)
       #print(x_sorted2)
       #print(y_sorted2)
       #print(z_sorted2)
       
       return x_sorted2,y_sorted2,z_sorted2

    
  
    #def rotoTraslate (self, x,y,x0,y0,phi):
    def rotoTraslate(self, x,y):
          
        x0=self.x0
        y0=self.y0
        phi=self.phi 
        dx = (x - x0)
        dy = (y - y0)
        xp = numpy.cos(phi)*dx + numpy.sin(phi)*dy
        yp = -numpy.sin(phi)*dx + numpy.cos(phi)*dy
        #aaa=[xp,yp]

        #ordinare per x crescente
        #xs,ys=self.ordinateByX(xp,yp,z)
        
                
        
        return xp,yp

    
