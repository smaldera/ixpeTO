from ROOT import *

import numpy as np
import math
import sys

def fit_histo_phi(h):


   print "inizio fit func "
   
   fitFunc2 = TF1("fitFunc2", "[0]+[1]*cos(x-[2])*cos(x-[2])", -TMath.Pi(), TMath.Pi()) #LEGGE di MALUS
      
   fitFunc2.SetParameter(0,1000)
   fitFunc2.SetParameter(1,1500)
   
   fitFunc2.SetParLimits(0,0,100000000)     # offset >0
   fitFunc2.SetParLimits(1,0,100000000)     # se c'e' bisogno di un'inversione, che la becchi con la fase [2]!!!
   print "prima del fit... "

   h.Fit(fitFunc2,"MR")
   
   print "fit fatto... "
   gStyle.SetOptStat("nem")
   gStyle.SetOptFit(100)
   gStyle.SetStatW(0.1)
   gStyle.SetStatH(0.09)
   
   c=TCanvas("c","",0)  # OKKIO!!! without this line it crashes at runtime!!!! ( why???)
   h.Draw("E1") #"E1" to show error bars
      
   reducedChiSquare = fitFunc2.GetChisquare()/fitFunc2.GetNDF() #getProb
   fitter = TVirtualFitter.GetFitter()
   covMatrix = fitter.GetCovarianceMatrix()
   print 'cov(A,B) = cov(0,1) = ', fitter.GetCovarianceMatrixElement(0,1)
   covAB = fitter.GetCovarianceMatrixElement(0,1)


   A = fitFunc2.GetParameter(0)
   B = fitFunc2.GetParameter(1)
   sA = fitFunc2.GetParError(0)
   sB = fitFunc2.GetParError(1)
   modulationFactor = fitFunc2.GetParameter(1)/(fitFunc2.GetParameter(1)+2*fitFunc2.GetParameter(0))
   sModulationFactor = math.sqrt( ((4*B*B*sA*sA)+(4*A*A*sB*sB)-(8*A*B*covAB) ) / ((2*A+B)*(2*A+B)*(2*A+B)*(2*A+B)) ) #error with covariance
     

   mu= [modulationFactor,sModulationFactor]
      
   return mu
   



if __name__ == '__main__':


  
#if number of args not correct
  if (len(sys.argv) != 3):
    print "    usage: python plot_modFactor.py   configFile  OutfileName.root"
    sys.exit()


print "reading root file list from: ",sys.argv[1]    


    
# read config file

miofile = open(sys.argv[1],'r') 
data = miofile.readlines()  # -> questo va bene per file piccoli!!! carica tutto il file sulla ram
n_points = len(data)
print data

print "n_points = ",n_points

mod_fact=np.array( [0.]*n_points,"f")
mod_factErr=np.array( [0.]*n_points,"f")
energy_Err=np.array( [0.]*n_points,"f")
energy=np.array( [0.]*n_points,"f")
mod_factDiff=np.array( [0.]*n_points,"f")
mod_factDiffErr=np.array( [0.]*n_points,"f")




for i in range (0, len(data)):  # for sulle linee del file

   line=data[i]
   print "processing file = ",line.split()[0]," energy[KeV]=  ",line.split()[1]
   

   filename= line.split()[0]
   energy[i]=float(line.split()[1])
   
   
   print "\n\n infile = ",filename
   f=TFile(filename)
   h=f.Get("5,5")
   modulationFactor=-1
   sModulationFactor=-1
   mu= fit_histo_phi(h)
   mod_fact[i]= mu[0]
   mod_factErr[i]=mu[1]
   energy_Err[i]=0.
    
   print "modFactor standard thersholds(5,5) = ",modulationFactor," +/- ",sModulationFactor
 
   xp=0
   yp=0
   
   try:
      cc=f.Get("cc")
   except:
      print "file = ",filename," no cc ... stop here \n\n"
      break

   try:
      maxP=cc.GetPrimitive("TMarker")
   except:
      print "file = ",filename," no TMarker... stiop here \n\n "
      break 

  
   xp=int (maxP.GetX())
   yp=int (maxP.GetY())
      
   nomeHist=str(xp)+","+str(yp)
   print "get histogrem =",nomeHist
   h=f.Get(nomeHist)
   muMax= fit_histo_phi(h)
   print "modFactor Max (",xp,",","yp",") = ",muMax[0]," +/- ",muMax[1]
 
   mod_factDiff[i]=(muMax[0] - mu[0] )*100.
   mod_factDiffErr[i]= math.sqrt( (muMax[1]*muMax[1])+(mu[1]*mu[1]) )*100.
      
   f.Close()
  
      
g= TGraphErrors(len(mod_fact),energy,mod_fact, energy_Err,mod_factErr)
g.SetName("mu_std")
g.SetTitle("mu_std (5,5)" )
g.GetXaxis().SetTitle("E [KeV]" )
g.GetYaxis().SetTitle("#mu" )


g.SetMarkerStyle(24)
c=TCanvas("c","",0)
g.Draw("ap")


gDiff= TGraphErrors(len(mod_fact),energy,mod_factDiff, energy_Err,mod_factDiffErr)
gDiff.SetName("mu_diff")
gDiff.SetMarkerStyle(20)
gDiff.SetTitle("mu_diff" )
gDiff.GetXaxis().SetTitle("E [KeV]" )
gDiff.GetYaxis().SetTitle("( #mu_{max}-#mu_{std} ) x 100" )


c2=TCanvas("c2","",0)
gDiff.Draw("ap")



print "\n\nopening out root file...",sys.argv[2]
fout=TFile(sys.argv[2],"recreate")
print "std mod factor plot saved  as mu_std"
g.Write()
print "mod factor difference  plot saved  as mu_diff"
gDiff.Write()
print "... all done!! "
fout.Close()
