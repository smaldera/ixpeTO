# classe per raggruppare le varie funzioni usate per fittare la distibuzione 1d di carica (ottenuta dalla proiezione lungo la traiettoria fittata...)
#

import ROOT



class fitDistrib:

   def  super_Fit2Gaussiane(self,h1):

       convp_original=self.myPeakFind(h1)
       redChi2_min=2e13
       convp_min=0
       for i in range (0,10):
              convp=convp_original-0.05+0.01*i
              #convp=convp_original
              f=self.Fit2Gaussiane_convp(h1,convp)

              if (f.GetNDF()>0):
                         redChi2=f.GetChisquare()/f.GetNDF()
                         
              else:
                        redChi2=1e12
              print("i=",i, " convp=",convp," (",convp_original,") redChi2=",redChi2)
              if redChi2<redChi2_min:
                     redChi2_min=redChi2
                     convp_min=convp

       print ("convp min = ",convp_min)
       return self.Fit2Gaussiane_convp(h1,convp_min)

                     
       


   def  Fit2Gaussiane_convp (self,h1,convp):

       G0 = ROOT.TF1 ("G0","gaus",-1,1)
       h1.Fit("G0","MEWWRN")
       mean0=G0.GetParameter(1)
       sigma0=G0.GetParameter(2)
       maxH = h1.GetBinCenter( h1.GetMaximumBin())
       #convp=-0.2
       print("convp = ",convp)

       G2 = ROOT.TF1 ("G2","gaus",maxH-1.0*sigma0,maxH+1.0*sigma0)
       G2.SetParameter(0,  h1.GetMaximumBin() )
       G2.SetParameter(1,maxH )
       G2.SetParameter(2, sigma0 )
       G2.SetParLimits(1, mean0-sigma0/2., mean0+sigma0/2.)
       h1.Fit("G2","MERwwNq")
         

       G1 = ROOT.TF1 ("G1","gaus",convp-0.05,convp+0.05)
       G1.SetParameter(0,3.33294e+02 )
       G1.SetParameter(1, convp )
       G1.SetParameter(2, 5.30043e-02  )       
       h1.Fit("G1","MERwwNq")

       
       #input("continue2?")
       
       Gsum=ROOT.TF1("Gsum","gaus(0)+gaus(3)",-0.4,0.4)
       #Gsum=ROOT.TF1("Gsum","gaus(0)+gaus(3)",-0.3,0.3)
       Gsum.SetParameter(0, G1.GetParameter(0))
       Gsum.SetParameter(1, G1.GetParameter(1))
       Gsum.SetParameter(2, G1.GetParameter(2))
       Gsum.SetParameter(3, G2.GetParameter(0))
       Gsum.SetParameter(4, G2.GetParameter(1))
       Gsum.SetParameter(5, G2.GetParameter(2))

       
       #Gsum.SetParLimits(0,  G1.GetParameter(0)-G1.GetParameter(0)*0.2,  G1.GetParameter(0)+ G1.GetParameter(0)*0.2)
       Gsum.SetParLimits(0,  0,  1000)
       #Gsum.SetParLimits(1,  G1.GetParameter(1)-G1.GetParameter(1)*0.8,  min(G1.GetParameter(1)+ G1.GetParameter(1)*0.8,0))
       Gsum.SetParLimits(1,  convp-0.03, convp+0.03 )
       #Gsum.SetParLimits(2,  G1.GetParameter(2)-G1.GetParameter(2)*0.5,  G1.GetParameter(2)+ G1.GetParameter(2)*0.5)

       Gsum.SetParLimits(3,  0,  1000)
       #Gsum.SetParLimits(3,  G2.GetParameter(0)-G2.GetParameter(0)*0.2,  G2.GetParameter(0)+ G2.GetParameter(0)*0.2)
       #Gsum.SetParLimits(4,  G2.GetParameter(1)-G2.GetParameter(1)*0.5,  G2.GetParameter(1)+ G2.GetParameter(1)*0.5)
       #Gsum.SetParLimits(5,  G2.GetParameter(2)-G2.GetParameter(2)*0.5,  G2.GetParameter(2)+ G2.GetParameter(2)*0.5)
       
       #h1.Fit("Gsum","LMERwwN")
       h1.Fit("Gsum","MERWNq")
       #Gsum.Draw("samel")

       if (Gsum.GetNDF()>0):
           redChi2=Gsum.GetChisquare()/Gsum.GetNDF()
           print("redChi2= ",redChi2)
       else:
           redChi2=1e14
         

       
       #Gsum=0 
       return Gsum


   def fitGaus_cutoff_old (self, h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print ("fit con cut-off")
        
        max = h.GetBinCenter( h.GetMaximumBin())
          
        convp=self.myPeakFind(h)
        G0 = ROOT.TF1 ("G0","gaus",max-0.1,max+0.1)
        G0.SetParameter(1,max)
        h.Fit("G0","LMERwwN")
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        #G0.Draw("samel")
        
                
        G1=ROOT.TF1 ("G1","gaus",convp-0.1,convp+0.1)
        G1.SetLineColor(4)
        h.Fit("G1","LMERwwN")
        #G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)
                
        f1= ROOT.TF1 ("f1","  ([0]*x-[1])*(1./  (1+exp( (x-[2])/[3])  ))", mean0-sigma0,max+0.2)
        f1.SetParameters(500,88,mean0,0.01)
        
        f1.SetParLimits(2, mean0-sigma0/2., mean0+sigma0)
        h.Fit("f1","LMERwwN")
       
        fsum=ROOT.TF1("fsum","gaus(0)+ (  ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))*(  ( ( ([3]*x-[4])*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  )))) >0  )  )",convp-0.2,max+0.2)
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)

        fsum.SetParLimits(0, 30,300)
        fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0) ) 
        fsum.SetParLimits(2, 0.02, sigma0)
        fsum.SetParLimits(5,  mean0-sigma0, mean0+sigma0)
        
        fsum.SetLineColor(2)
        fsum.SetLineWidth(4)
        h.Fit("fsum","LMERwwN")

        return fsum




   def fitExp_cutoff_old (self, h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print ("fit con gaus + exp+ cut-off")
        max = h.GetBinCenter( h.GetMaximumBin())
        convp=self.myPeakFind(h)
        print ("distConv= ",self.distConv)


        G0 = ROOT.TF1 ("G0","gaus",max-0.1,max+0.1)
        h.Fit("G0","LMERN")
        G0.SetParameter(1,max)
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        #G0.Draw("samel")

        
        
        #OK !!
        G1=ROOT.TF1 ("G1","gaus",convp-0.1,convp+0.1)
        G1.SetLineColor(4)
        h.Fit("G1","LMERN")
        #G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)

        f1= ROOT.TF1 ("f1","  (exp([0]+[1]*x))*(1./  (1+exp( (x-[2])/[3])  ))", convp+0.1,max+0.2)
        f1.SetParameters(2.3,4,max,0.01)

        f1.SetParLimits(2, max-sigma0/2., max+sigma0)
        #f1.FixParameter(2, mean0)
        h.Fit("f1","LMERN")
        f1.SetLineColor(3)
        #f1.Draw("samel")
        
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

        
        fsum.SetParLimits(0, 15,200)
       
        fsum.SetParLimits(1, convp-0.1 , min(convp+0.2,max-sigma0) ) 
        fsum.SetParLimits(2, 0.02, sigma0) #!!!!!!!!!!!!!
        fsum.SetParLimits(5,  mean0-sigma0, mean0+sigma0)
        
        fsum.SetLineColor(2)
        fsum.SetLineWidth(4)
        h.Fit("fsum","LMERN")

        #if (fsum.GetNDF()>0):
        #    self.redChi2=fsum.GetChisquare()/fsum.GetNDF()
        #else:
        #    self.redChi2=1000 
            
        return fsum


     
   def fit2Gaussiane (self,h1):

       print("fit con doppia gaussiana")
       G0 = ROOT.TF1 ("G0","gaus",-1,1)
       h1.Fit("G0","MEWWRN")
       mean0=G0.GetParameter(1)
       sigma0=G0.GetParameter(2)
       maxH = h1.GetBinCenter( h1.GetMaximumBin())
       convp=self.myPeakFind(h1)

       G2 = ROOT.TF1 ("G2","gaus",maxH-1.5*sigma0,maxH+1.5*sigma0)
       G2.SetParameter(0,  h1.GetMaximumBin() )
       G2.SetParameter(1,maxH )
       G2.SetParameter(2, sigma0 )

       G2.SetParLimits(1, mean0-sigma0/2., mean0+sigma0/2.)
       #G2.FixParameter(1)
       #G2.FixParameter(2)
       h1.Fit("G2","MERwwN")
         

       
       G1 = ROOT.TF1 ("G1","gaus",convp-0.05,convp+0.05)
       G1.SetParameter(0,3.33294e+02 )
       G1.SetParameter(1, convp )
       G1.SetParameter(2, 5.30043e-02  )       
       h1.Fit("G1","MERwwN")
       
       
       Gsum=ROOT.TF1("Gsum","gaus(0)+gaus(3)",-0.4,0.4)
       Gsum.SetParameter(0, G1.GetParameter(0))
       Gsum.SetParameter(1, G1.GetParameter(1))
       Gsum.SetParameter(2, G1.GetParameter(2))
       Gsum.SetParameter(3, G2.GetParameter(0))
       Gsum.SetParameter(4, G2.GetParameter(1))
       Gsum.SetParameter(5, G2.GetParameter(2))
       
       #Gsum.SetParLimits(0,  G1.GetParameter(0)-G1.GetParameter(0)*0.2,  G1.GetParameter(0)+ G1.GetParameter(0)*0.2)
       Gsum.SetParLimits(0,  0,  1000)
       Gsum.SetParLimits(1,  G1.GetParameter(1)-G1.GetParameter(1)*0.8,  G1.GetParameter(1)+ G1.GetParameter(1)*0.8)
       #Gsum.SetParLimits(2,  G1.GetParameter(2)-G1.GetParameter(2)*0.5,  G1.GetParameter(2)+ G1.GetParameter(2)*0.5)

       Gsum.SetParLimits(3,  0,  1000)
       #Gsum.SetParLimits(3,  G2.GetParameter(0)-G2.GetParameter(0)*0.2,  G2.GetParameter(0)+ G2.GetParameter(0)*0.2)
       #Gsum.SetParLimits(4,  G2.GetParameter(1)-G2.GetParameter(1)*0.5,  G2.GetParameter(1)+ G2.GetParameter(1)*0.5)
       #Gsum.SetParLimits(5,  G2.GetParameter(2)-G2.GetParameter(2)*0.5,  G2.GetParameter(2)+ G2.GetParameter(2)*0.5)
       
       
       h1.Fit("Gsum","LMERwwN")
       #Gsum.Draw("samel")
       if (Gsum.GetNDF()>0):
           redChi2=Gsum.GetChisquare()/Gsum.GetNDF()
           print("redChi2= ",redChi2)
       else:
           redChi2=1000 
              
       return Gsum
       
############################3
   def fitExp_cutoff (self,h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print ("fit con doppio  exp+ cut-off")
        maxH = h.GetBinCenter( h.GetMaximumBin())
        convp=myPeakFind(h)
        #convp=-0.12
       


        G0 = ROOT.TF1 ("G0","gaus",maxH-0.5,+maxH+0.5)
        
        G0.SetParameter(1,maxH)
        h.Fit("G0","LMERWWN")
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        #G0.Draw("samel")
        G1=ROOT.TF1 ("G1","(exp([0]+[1]*x))*(1./  (1+exp( (x-[2])/[3])  ))",convp-0.10,convp+0.05)
        G1.SetLineColor(4)
        G1.SetParameter(0, 16 )
        G1.SetParameter(1, -10 )
        
        G1.SetParameter(2, convp )
        G1.SetParameter(3, 0.5 )
        G1.SetParLimits(2, convp-0.03,convp+0.03 ) 
      
        h.Fit("G1","LMERwN")
        #G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)
        #G1.Draw("samel")
        
        
        f1= ROOT.TF1 ("f1","  (exp([0]+[1]*x) )*(1./  (1+exp( (x-[2])/[3])  ))", -0.5, 0.5  )
        f1.SetParameters(2.3,4, maxH-0.01,0.01)
        
      
        h.Fit("f1","LMRWWN")
        f1.SetLineColor(3)
        #f1.Draw("samel")

        fsum=ROOT.TF1("fsum","  (exp([0]+[1]*x))*(1./  (1+exp( (x-[2])/[3])  ))   + (  ( exp([4]+[5]*x)*(1-exp(-x/[8]) )*(1./  (1+exp( (x-[6])/[7])  ))))",-0.4,0.4)
      
        #ok
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        fsum.SetParameter(3, G1.GetParameter(3))
        
        fsum.SetParameter(4, f1.GetParameter(0))
        fsum.SetParameter(5, f1.GetParameter(1))
        fsum.SetParameter(6, f1.GetParameter(2))
        fsum.SetParameter(7, f1.GetParameter(3))
        fsum.SetParameter(8, 0.1)


        fsum.SetParLimits(0, G1.GetParameter(0)- G1.GetParameter(0)*0.2,  G1.GetParameter(0)+ G1.GetParameter(0)*0.2 )
        #fsum.SetParLimits(0, 15,200)
        fsum.SetParLimits(1,  G1.GetParameter(1)-G1.GetParameter(1)*0.2,  G1.GetParameter(1)+ G1.GetParameter(1)*0.2)
        fsum.SetParLimits(2,  G1.GetParameter(2)-G1.GetParameter(2)*0.2,  G1.GetParameter(2)+ G1.GetParameter(2)*0.2)
        fsum.SetParLimits(3,  G1.GetParameter(3)-G1.GetParameter(3)*0.2,  G1.GetParameter(3)+ G1.GetParameter(3)*0.21)



        h.Fit("fsum","LMERWWN")
        if (fsum.GetNDF()>0):
           redChi2=fsum.GetChisquare()/fsum.GetNDF()
           print("redChi2= ",redChi2)
        else:
           redChi2=1000 
       
        return fsum


   
   def fitExpCutoff_gaus (self,h):
        ###         
        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        ###
          
        print ("fit con gaus + exp+ cut-off")
        maxH = h.GetBinCenter( h.GetMaximumBin())
        convp=self.myPeakFind(h)

        G0 = ROOT.TF1 ("G0","gaus",maxH-0.5,+maxH+0.5)
        G0.SetParameter(1,maxH)
        h.Fit("G0","LMERWWN")
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        #G0.Draw("samel")

                
        G1=ROOT.TF1 ("G1","gaus",convp-0.05,convp+0.05)
        G1.SetLineColor(4)
        G1.SetParameter(2, convp )
        G1.SetParameter(3, 0.5 )
        G1.SetParLimits(2, convp-0.03,convp+0.03 ) 
     
        h.Fit("G1","LMERwN")
        #G1.Draw("samel")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)
        #G1.Draw("samel")
       
        
        f1= ROOT.TF1 ("f1","  (exp([0]+[1]*x) )*(1./  (1+exp( (x-[2])/[3])  ))", -0.5, 0.5  )
        f1.SetParameters(2.3,4, maxH-0.01,0.01)

        h.Fit("f1","LMRWWN")
        #f1.SetLineColor(3)
        #f1.Draw("samel")
        
        fsum=ROOT.TF1("fsum"," gaus(0)+ (  ( exp([3]+[4]*x)*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  ))))",-0.4,0.4)
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        
        
        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)

        #ok
        #fsum.SetParLimits(0, G1.GetParameter(0)- G1.GetParameter(0)*0.2,  G1.GetParameter(0)+ G1.GetParameter(0)*0.2 )
        #fsum.SetParLimits(0, 1,1000)
        fsum.SetParLimits(1,  G1.GetParameter(1)-G1.GetParameter(1)*0.5,  G1.GetParameter(1)+ G1.GetParameter(1)*0.5)
        #fsum.SetParLimits(1,  convp-0.1, convp+0.1 )
       
        #fsum.SetParLimits(2,  G1.GetParameter(2)-G1.GetParameter(2)*0.2,  G1.GetParameter(2)+ G1.GetParameter(2)*0.2)               
        #fsum.SetParLimits(5,  f1.GetParameter(2)-f1.GetParameter(2)*0.9,  f1.GetParameter(2)+ f1.GetParameter(2)*0.9)
        
        
        h.Fit("fsum","LMERWWN")
        if (fsum.GetNDF()>0):
           redChi2=fsum.GetChisquare()/fsum.GetNDF()
           print("redChi2= ",redChi2)
        else:
           redChi2=1000 
         
        return fsum



############################3
   def fitExpCutoff_landau (self,h):

        #ottimizzato per dividiBins=1 ,filtro: M=18 cutoff=7e7
        
        print ("fit con gaus + exp+ cut-off")
        maxH = h.GetBinCenter( h.GetMaximumBin())
        convp=self.myPeakFind(h)

        G0 = ROOT.TF1 ("G0","gaus",maxH-0.5,+maxH+0.5)
        G0.SetParameter(1,maxH)
        h.Fit("G0","LMERWWN")
        mean0=G0.GetParameter(1)
        sigma0=G0.GetParameter(2)
        #G0.Draw("samel")

        
        
        #OK !!
        G1=ROOT.TF1 ("G1","landau",convp-0.05,convp+0.05)
        G1.SetLineColor(4)
        
        G1.SetParameter(2, convp )
        G1.SetParameter(3, 0.5 )
        G1.SetParLimits(2, convp-0.03,convp+0.03 ) 
       
        h.Fit("G1","LMERwN")
        mean1=G1.GetParameter(1)
        sigma1=G1.GetParameter(2)
        #G1.Draw("samel")
        
        f1= ROOT.TF1 ("f1","  (exp([0]+[1]*x) )*(1./  (1+exp( (x-[2])/[3])  ))", -0.5, 0.5  )
        f1.SetParameters(2.3,4, maxH-0.01,0.01)
        h.Fit("f1","LMRWWN")
        #f1.SetLineColor(3)
        #f1.Draw("samel")

        
        
        fsum=ROOT.TF1("fsum"," landau(0)+ (  ( exp([3]+[4]*x)*(1-exp(-x/[7]) )*(1./  (1+exp( (x-[5])/[6])  ))))",-0.4,0.4)
        fsum.SetParameter(0, G1.GetParameter(0))
        fsum.SetParameter(1, G1.GetParameter(1))
        fsum.SetParameter(2, G1.GetParameter(2))
        fsum.SetParameter(3, f1.GetParameter(0))
        fsum.SetParameter(4, f1.GetParameter(1))
        fsum.SetParameter(5, f1.GetParameter(2))
        fsum.SetParameter(6, f1.GetParameter(3))
        fsum.SetParameter(7, 0.1)
        
        #fsum.SetLineWidth(4)
        h.Fit("fsum","LMERWWN")
        if (fsum.GetNDF()>0):
           redChi2=fsum.GetChisquare()/fsum.GetNDF()
           print("redChi2= ",redChi2)
        else:
           redChi2=1000 
       
        return fsum

 ######################################

   def myPeakFind(self,h):
  
    # minPeakIId=-2.5
     minPeakIId=-1.8 
     
     nBins=h.GetNbinsX()
     hDiff=h.Clone()
     for i in range (2,nBins):
           counts_i=h.GetBinContent(i)
           counts_prev=h.GetBinContent(i-1)
           diff=counts_i-counts_prev
          # print ("i=",i, " diff = ",diff)
           hDiff.SetBinContent(i,diff)
           
     xZero=-1e10   
     for i in range (2,nBins):      
          counts_i=hDiff.GetBinContent(i)
          counts_prev=hDiff.GetBinContent(i-1)    
          if ( counts_i<7 and  counts_prev>7 and  h.GetBinContent(i)>0    ):
                 xZero=(h.GetBinCenter(i)+h.GetBinCenter(i-1))/2.
                 print("xero found... ",xZero)
                 break

     hDiff2=h.Clone()
     for i in range (2,nBins):
           counts_i=hDiff.GetBinContent(i)
           counts_prev=hDiff.GetBinContent(i-1)
           diff2=counts_i-counts_prev
          # print ("i=",i, " diff = ",diff)
           hDiff2.SetBinContent(i,diff2)       

     xZero2=-1e10   
     for i in range (2,nBins):
          
          xZero2=-1e10
          counts_i=hDiff2.GetBinContent(i)
          counts_prev=hDiff2.GetBinContent(i-1)    
          #if ( counts_i<0 and  counts_prev>0 and  h.GetBinContent(i)>0    ):
          if ( counts_i<0 and  counts_prev>0 and  h.GetBinContent(i)>50    ):

             print("trovato zero2 a i=",i, " ",h.GetBinCenter(i) )
             # cerco minimo successivo
             primoMinimo=0
             for jj in range (i,i+30):
                  #primoMinimo=0
                  print(" -- cerco minimo successivo, jj = ",jj," ",h.GetBinCenter(jj))
                  if (hDiff2.GetBinContent(jj)- hDiff2.GetBinContent(jj-1 )<-0.05) and (hDiff2.GetBinContent(jj)- hDiff2.GetBinContent(jj+1 )<-0.05):
                 # if (hDiff2.GetBinContent(jj)- hDiff2.GetBinContent(jj-1 )<-0.000005) and (hDiff2.GetBinContent(jj)- hDiff2.GetBinContent(jj+1 )<-0.00005):
                     
                     primoMinimo=hDiff2.GetBinContent(jj)
                     print("    ===>>>  TROVATO --min successivo, jj = ",jj," min = ",hDiff2.GetBinContent(jj), " ",h.GetBinCenter(jj) )
                     break
                  
                     
             if primoMinimo<minPeakIId or ( primoMinimo==0 and hDiff2.GetBinContent(i+30-1)<minPeakIId ):   
                  xZero2=(h.GetBinCenter(i)+h.GetBinCenter(i-1))/2.
                  print("xero 2found... ",xZero2)
                  break       


     #c2=ROOT.TCanvas("c2","",0)       
     #hDiff.Draw()
     #hDiff2.SetLineColor(2)
     #hDiff2.Draw("same")
     #input("contiunue?")
     

            
     return xZero2




####################################################

def testall():
   files=["","2","3","4","6","7","8","9","11","12","13","14","15","16"]
   #files=["13","14","15","16"]
   #files=["17"]
   for num in files:
      filename="pippo"+num+".root"
      print ("")
      print ("=======================================")
      print ("processing file:",filename)
      print ("=======================================")
      print ("") 
      inRootFile=ROOT.TFile(filename,"open")
      h=inRootFile.Get("h1L")
      h2=h.Clone()
      h2.Draw()

      a=fitDistrib()
      a.super_Fit2Gaussiane(h2) # abbastanza buono per E=5KeV

      input("continue3?")

      
if __name__ == '__main__':


   testall()
   """
   inRootFile=ROOT.TFile("pippo.root","open")
   h=inRootFile.Get("h1L")
   h2=h.Clone()
   h2.Draw()

   a=fitDistrib()
   a.super_Fit2Gaussiane(h2)
   #fit2Gaussiane(h2) # ~ ok
   """

   #fitExp_cutoff(h2) # circa ok
   #fitExpCutoff_gaus(h2)


   #fitExpCutoff_landau(h2) noooo


   #fitLin_cutoff(h2) # ELIMITANTO

   #myPeakFind(h2)
