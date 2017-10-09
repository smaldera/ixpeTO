import ROOT
import math

def fit_modFact(h1):
    
    
    fitFunc = ROOT.TF1("fitFunc", "[0]+[1]*cos(x-[2])*cos(x-[2])", -ROOT.TMath.Pi(), ROOT.TMath.Pi()) #LEGGE di MALUS
    fitFunc.SetParLimits(0,0,100000000)     # offset >0
    fitFunc.SetParLimits(1,0,100000000)     # se c'e' bisogno di un'inversione, che la becchi con la fase [2]!!! 
    h1.Fit("fitFunc","MR")

    reducedChiSquare = fitFunc.GetChisquare()/fitFunc.GetNDF() #getProb
    probability = fitFunc.GetProb()

    #get covariance matrix
    fitter = ROOT.TVirtualFitter.GetFitter()
    covMatrix = fitter.GetCovarianceMatrix()
    
    print 'cov(A,B) = cov(0,1) = ', fitter.GetCovarianceMatrixElement(0,1)
    covAB = fitter.GetCovarianceMatrixElement(0,1)
    
    A = fitFunc.GetParameter(0)
    B = fitFunc.GetParameter(1)
    sA = fitFunc.GetParError(0)
    sB = fitFunc.GetParError(1)
    modulationFactor = B/(B+2*A)
    sModulationFactor = math.sqrt( ((4*B*B*sA*sA)+(4*A*A*sB*sB)-(8*A*B*covAB) ) / ((2*A+B)*(2*A+B)*(2*A+B)*(2*A+B)) ) #error with covariance CONTROLLARE FORMULA!!!!!

    mu=[modulationFactor,sModulationFactor]
    return mu




def dist90 (hist):
    n=hist.GetEntries()
    nBins=hist.GetNbinsX()
    d90=-1.
    d68=-1.
    
    for i in range (0,nBins):
           ratio=hist.Integral(0,i)/float(n)
           #print "i =", i," ratio=",ratio, " enttries = ",n
           if (ratio>=0.68):
                d68=hist.GetBinCenter(i)
               
           if (ratio>=0.9):
                 d90=hist.GetBinCenter(i)
                 #print " d90 = ",d90
                 break


    dists=[d90, d68]         
             
    return dists


 

def readScan(inputFile, outFile):
    inRootFile= ROOT.TFile(inputFile,"open")


    h2d=ROOT.TH2F("h2","dMean",16, 0,8,16,0,8)
    h2d_offset=ROOT.TH2F("h2_offset","mean_offset",16, 0,8,16,0,8)
    h2d90=ROOT.TH2F("h2","d90",16, 0,8,16,0,8)
    h2d68=ROOT.TH2F("h2","d68",16, 0,8,16,0,8)
    


    h2d.GetXaxis().SetTitle("Rmin")
    h2d.GetYaxis().SetTitle("Rmax")
    h2d90.GetXaxis().SetTitle("Rmin")
    h2d90.GetYaxis().SetTitle("Rmax")
    h2d68.GetXaxis().SetTitle("Rmin")
    h2d68.GetYaxis().SetTitle("Rmax")
    h2d_offset.GetXaxis().SetTitle("Rmin")
    h2d_offset.GetYaxis().SetTitle("Rmax")
    
    tree=inRootFile.Get("tree_rec")

    print "N entries = ",tree.GetEntries()

    hDist=ROOT.TH1F("hDist","MC - Rec conversion point dist.",1000,0,1)
    hDiffXrt=ROOT.TH1F("hDiffXrt","",2000,-1,1)

    minR_std=1.5
    maxR_std=3.5

    minr_final=-1
    maxr_final=-1
    minr_off=-1
    maxr_off=-1
    minr_d90=-1
    maxr_d90=-1
    minr_d68=-1
    maxr_d68=-1
   

    
    minOffset=10000
    minD=100000
    minD90=100000
    minD68=100000
   
    
    stdD=-1
    stdOffset=-1
    stdD90=-1
    stdD68=-1
    
    
    
    for kk in range(0,10):
        raggioMin=float(0.5*kk)
        for jj in range (0,14):
            raggioMax=raggioMin+float(0.5*jj)
            if raggioMax> 8:
                break


            cut='raggioMin=='+str(raggioMin)+'&& raggioMax=='+str(raggioMax)
            print "cut = ",cut
            tree.Draw("dist2D>>hDist",cut)
            d90=dist90(hDist)[0]
            d68=dist90(hDist)[1]
            dMean=hDist.GetMean()
            
            tree.Draw("xRec_rt>>hDiffXrt",cut)
            #offset=((hDiffXrt.GetMean()**2)+(hDiffXrt.GetRMS()**2))**0.5
            offset=hDiffXrt.GetMean()
           
           
            h2d.Fill(raggioMin,raggioMax,dMean)
            h2d90.Fill(raggioMin,raggioMax,d90)
            h2d68.Fill(raggioMin,raggioMax,d68)
            h2d_offset.Fill(raggioMin,raggioMax,offset)
            

            if (raggioMin== minR_std and raggioMax== maxR_std):
                 stdD=dMean
                 stdOffset=offset
                 stdD90=d90
                 stdD68=d68
            
           
            if offset<minOffset and hDiffXrt.GetEntries>100:
                minr_off=raggioMin
                maxr_off=raggioMax
                minOffset= offset

            if  dMean<minD and hDist.GetEntries>100:
                minr_final=raggioMin
                maxr_final=raggioMax
                minD=dMean 

                
                
            if d90<minD90 and hDist.GetEntries>100:
                minr_d90=raggioMin
                maxr_d90=raggioMax
                minD90= d90

            if d68<minD68 and hDist.GetEntries>100:
                minr_d68=raggioMin
                maxr_d68=raggioMax
                minD68= d68

                
                
           
            hDiffXrt.Reset()
            hDist.Reset()
           #break
      #break
      
      


    
    outString= str(minr_final)+' '+ str(maxr_final)+' '+str(minD)+' '+str(minr_d90)+'  '+str(maxr_d90)+' '+str(minD90)+' '+str(minr_d68)+'  '+str(maxr_d68)+' '+str(minD68)+' '+str(minr_off)+' '+ str(maxr_off)+' '+ str(minOffset)+' ' +str(stdD)+'  '+str(stdD90)+'  '+str(stdD68)+' '+ str(stdOffset)  

      
    miofile = open('out.txt','w')       # apre il file in scrittura
    miofile.write(outString)

 

      

      ##############


      
    hDist_std=ROOT.TH1F("hDist_std","MC - rec conversion point dist.",1000,0,1)
    hDist_max=ROOT.TH1F("hDist_max","MC - rec conversion point dist.",1000,0,1)
    hDist90=ROOT.TH1F("hDist90","MC - rec conversion point dist.",1000,0,1)
    hDist68=ROOT.TH1F("hDist68","MC - rec conversion point dist.",1000,0,1)
    hDistOffset=ROOT.TH1F("hDistOffset","MC - rec conversion point dist.",1000,0,1)
   
    

    hPhi=ROOT.TH1F("hPhi","", 90, -3.14, 3.14)
    hPhi_std=ROOT.TH1F("hPhi_std","", 90, -3.14, 3.14)
    hPhi_offset=ROOT.TH1F("hPhi_offset","", 90, -3.14, 3.14)
    
    hPhiDiff=ROOT.TH1F("hPhiDiff","", 360, -3.14, 3.14)
    hPhiDiff_std=ROOT.TH1F("hPhiDiff_std","", 360, -3.14, 3.14)
    
    hDist_max.GetXaxis().SetTitle("2D_dist")
    hDist_std.GetXaxis().SetTitle("2D_dist")
    hPhi_std.GetXaxis().SetTitle("#phi")
    hPhiDiff_std.GetXaxis().SetTitle("#Delta#phi")
    hPhiDiff.GetXaxis().SetTitle("#Delta#phi")
      
    hPhi_std.SetLineColor(4)
    hPhi.SetLineColor(2)
    hPhiDiff_std.SetLineColor(4)
    hPhiDiff.SetLineColor(2)
    hDist_std.SetLineColor(4)
    hDist_max.SetLineColor(2)
      


    cut='raggioMin=='+str(minR_std)+'&& raggioMax=='+str(maxR_std)
    tree.Draw("dist2D>>hDist_std",cut)
    tree.Draw("phiRec>>hPhi_std",cut)
    tree.Draw("phiRec-phiMC>>hPhiDiff_std",cut)


    
    # massimo valore medio
    cut='raggioMin=='+str(minr_final)+'&& raggioMax=='+str(maxr_final)
    tree.Draw("dist2D>>hDist_max",cut)
    tree.Draw("phiRec-phiMC>>hPhiDiff",cut)
    tree.Draw("phiRec>>hPhi",cut)

    #massimo d90:
    cut='raggioMin=='+str(minr_d90)+'&& raggioMax=='+str(maxr_d90)
    tree.Draw("dist2D>>hDist90",cut)
    
    #massimo d68:
    cut='raggioMin=='+str(minr_d68)+'&& raggioMax=='+str(maxr_d68)
    tree.Draw("dist2D>>hDist68",cut)
    
    #massimo offset:
    cut='raggioMin=='+str(minr_off)+'&& raggioMax=='+str(maxr_off)
    tree.Draw("dist2D>>hDistOffset",cut)
    tree.Draw("phiRec>>hPhi_offset",cut)
    

    
    

    mu_std= fit_modFact(hPhi_std)
    mu_max= fit_modFact(hPhi)
    mu_maxOff= fit_modFact(hPhi_offset)

    print "Mu std = ",mu_std[0]," +- ",mu_std[1]
    print "Mu max = ",mu_max[0]," +- ",mu_max[1]
    print "Mu maxOff = ",mu_maxOff[0]," +- ",mu_maxOff[1]
    
    
    markerMax=ROOT.TMarker(minr_final+0.25,maxr_final+0.25,20)


    c1=ROOT.TCanvas("c1","",0)
    h2d.Draw("colz")
    markerMax.Draw("samep")  

    c2=ROOT.TCanvas("c2","",0)
    h2d_offset.Draw("colz")
   
    
    #===== canvas 3
    c3=ROOT.TCanvas("c3","",0)
    c3.Divide(2,2)

    c3.cd(1) #------------------------
    h2d.Draw('colz')
    markerMax.Draw("samep")  
 
    c3.cd(2) #---------------------------
    hDist_max.Draw()
    hDist_std.Draw("sames")
    hDist90.Draw("sames")
    hDist68.Draw("sames")
    hDistOffset.Draw("sames")
    
    
    leg=ROOT.TLegend(0.4,0.6,0.7,0.8)
    leg.AddEntry(hDist,"dMEan max","l")
    leg.AddEntry(hDist90,"d90 max","l")
    leg.AddEntry(hDist68,"d68 max","l")
    leg.AddEntry(hDistOffset,"dOff max","l")   
    leg.AddEntry(hDist_std,"d90 std","l")
    

    
    leg.Draw()
      
      
    c3.cd(3) #------------------------------

    hPhi_std.Draw()
    hPhi.Draw("sames")
    hPhi_offset.Draw("sames")
   
    leg2=ROOT.TLegend(0.,0.85,0.6,0.95)
    legend='#phi max, #mu='+str(mu_max[0])+' #pm '+str(mu_max[1])
    legend_std='#phi std, #mu='+str(mu_std[0])+' #pm '+str(mu_std[1])
    legend_off='#phi maxOff, #mu='+str(mu_maxOff[0])+' #pm '+str(mu_maxOff[1])
   
    
    leg2.AddEntry(hPhi,legend,"l")
    leg2.AddEntry(hPhi_std,legend_std,"l")
    leg2.AddEntry(hPhi_offset,legend_off,"l")
  
    leg2.Draw()
      
    c3.cd(4) #------------------------------

    hPhiDiff.Draw()
    hPhiDiff_std.Draw("sames")
    leg3=ROOT.TLegend(0.15,0.7,0.4,0.85)
    leg3.AddEntry(hPhiDiff,"#phi diff max","l")
    leg3.AddEntry(hPhiDiff_std,"#phi diff std ","l")
    leg3.Draw()
    
      


    outRootFile=ROOT.TFile(outFile,"recreate")

    hDist_max.Write()
    hDist_std.Write()
    hPhi.Write()
    hPhi_std.Write()
    hPhiDiff.Write()
    hPhiDiff_std.Write()
    markerMax.Write()
    hDist90.Write()
    hDist68.Write()
    hDistOffset.Write()

    h2d.Write()
    h2d90.Write()
    h2d68.Write()
    h2d_offset.Write()
    
    c1.Write()
    c2.Write()
    c3.Write()

    outRootFile.Close()

    print "min = ",minr_final, " max = ",maxr_final," minD= ",minD,' ( std = ',stdD,')'
    print  " minOFF = ",minr_off, " max = ",maxr_off," minOffset= ",minOffset,' ( std = ',stdOffset,')'
    print  " minD90 = ",minr_d90, " maxD90 = ",maxr_d90," minD90= ",minD90
    print  " minD68 = ",minr_d68, " maxD68 = ",maxr_d68," minD68= ",minD68
    
    


     

if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('inputFile', type=str,  help='the input root file')
    parser.add_argument('outFile', type=str,  help='the output root file')
    args = parser.parse_args()
    
    readScan(args.inputFile,args.outFile)
