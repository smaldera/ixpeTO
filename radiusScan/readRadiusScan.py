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
    
    for i in range (0,nBins):
           ratio=hist.Integral(0,i)/float(n)
           #print "i =", i," ratio=",ratio, " enttries = ",n
           if (ratio>=0.9):
                 d90=hist.GetBinCenter(i)
                 #print " d90 = ",d90
                 break
    return d90


 

def readScan(inputFile, outFile):
    inRootFile= ROOT.TFile(inputFile,"open")


    h2d=ROOT.TH2F("h2","d90",16, 0,8,16,0,8)
    h2d_offset=ROOT.TH2F("h2_offset","mean_offset",16, 0,8,16,0,8)
    h2d.GetXaxis().SetTitle("Rmin")
    h2d.GetYaxis().SetTitle("Rmax")
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

    minOffset=10000
    minD=100000

    stdD=-1
    stdOffset=-1

    
    for kk in range(0,10):
        raggioMin=float(0.5*kk)
        for jj in range (0,14):
            raggioMax=raggioMin+float(0.5*jj)
            if raggioMax> 8:
                break


            cut='raggioMin=='+str(raggioMin)+'&& raggioMax=='+str(raggioMax)
            print "cut = ",cut
            tree.Draw("dist2D>>hDist",cut)
            d90=dist90(hDist)
            tree.Draw("xRec_rt>>hDiffXrt",cut)
            offset=((hDiffXrt.GetMean()**2)+(hDiffXrt.GetRMS()**2))**0.5
           
            h2d.Fill(raggioMin,raggioMax,d90)
            h2d_offset.Fill(raggioMin,raggioMax,offset)

            if (raggioMin== minR_std and raggioMax== maxR_std):
                 stdD=d90
                 stdOffset=offset
                
            
           
            if offset<minOffset and hDiffXrt.GetEntries>100:
                minr_off=raggioMin
                maxr_off=raggioMax
                minOffset= offset
           
            if d90<minD and hDist.GetEntries>100:
                minr_final=raggioMin
                maxr_final=raggioMax
                minD= d90
           
            hDiffXrt.Reset()
            hDist.Reset()
           #break
      #break
      
      


    print "min = ",minr_final, " max = ",maxr_final," minD= ",minD, " minOFF = ",minr_off, " max = ",maxr_off," minOffset= ",minOffset, " stdD = ", stdD, " stdOffset =",  stdOffset  
    outString= str(minr_final)+' '+ str(maxr_final)+' '+str(minD)+' '+str(minr_off)+' '+ str(maxr_off)+' '+ str(minOffset)+' ' +str(stdD)+' ' + str(stdOffset)  

      
    miofile = open('out.txt','w')       # apre il file in scrittura
    miofile.write(outString)

 

      

      ##############

    
    hDist_std=ROOT.TH1F("hDist_std","MC - rec conversion point dist.",1000,0,1)
    hPhi=ROOT.TH1F("hPhi","", 90, -3.14, 3.14)
    hPhi_std=ROOT.TH1F("hPhi_std","", 90, -3.14, 3.14)
    hPhiDiff=ROOT.TH1F("hPhiDiff","", 360, -3.14, 3.14)
    hPhiDiff_std=ROOT.TH1F("hPhiDiff_std","", 360, -3.14, 3.14)
    
    hDist.GetXaxis().SetTitle("2D_dist")
    hDist_std.GetXaxis().SetTitle("2D_dist")
    hPhi_std.GetXaxis().SetTitle("#phi")
    hPhiDiff_std.GetXaxis().SetTitle("#Delta#phi")
    hPhiDiff.GetXaxis().SetTitle("#Delta#phi")
      
    hPhi_std.SetLineColor(4)
    hPhi.SetLineColor(2)
    hPhiDiff_std.SetLineColor(4)
    hPhiDiff.SetLineColor(2)
    hDist_std.SetLineColor(4)
    hDist.SetLineColor(2)
      


    cut='raggioMin=='+str(minR_std)+'&& raggioMax=='+str(maxR_std)
    tree.Draw("dist2D>>hDist_std",cut)
    tree.Draw("phiRec>>hPhi_std",cut)
    tree.Draw("phiRec-phiMC>>hPhiDiff_std",cut)


    cut='raggioMin=='+str(minr_final)+'&& raggioMax=='+str(maxr_final)
    tree.Draw("dist2D>>hDist",cut)
    tree.Draw("phiRec-phiMC>>hPhiDiff",cut)
    tree.Draw("phiRec>>hPhi",cut)


    mu_std= fit_modFact(hPhi_std)
    mu_max= fit_modFact(hPhi)

    print "Mu std = ",mu_std[0]," +- ",mu_std[1]
    print "Mu max = ",mu_max[0]," +- ",mu_max[1]

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
    hDist.Draw()
    hDist_std.Draw("sames")
    leg=ROOT.TLegend(0.4,0.6,0.7,0.8)
    leg.AddEntry(hDist,"d90 max","l")
    leg.AddEntry(hDist_std,"d90 std","l")
    leg.Draw()
      
      
    c3.cd(3) #------------------------------

    hPhi_std.Draw()
    hPhi.Draw("sames")
    leg2=ROOT.TLegend(0.,0.85,0.6,0.95)
    legend='#phi max, #mu='+str(mu_max[0])+' #pm '+str(mu_max[1])
    legend_std='#phi std, #mu='+str(mu_std[0])+' #pm '+str(mu_std[1])
    
    leg2.AddEntry(hPhi,legend,"l")
    leg2.AddEntry(hDist_std,legend_std,"l")
    leg2.Draw()
      
    c3.cd(4) #------------------------------

    hPhiDiff.Draw()
    hPhiDiff_std.Draw("sames")
    leg3=ROOT.TLegend(0.15,0.7,0.4,0.85)
    leg3.AddEntry(hPhiDiff,"#phi diff max","l")
    leg3.AddEntry(hPhiDiff_std,"#phi diff std ","l")
    leg3.Draw()
    
      


    outRootFile=ROOT.TFile(outFile,"recreate")

    hDist.Write()
    hDist_std.Write()
    hPhi.Write()
    hPhi_std.Write()
    hPhiDiff.Write()
    hPhiDiff_std.Write()
    markerMax.Write()
    
    c1.Write()
    c2.Write()
    c3.Write()

    outRootFile.Close()



if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('inputFile', type=str,  help='the input root file')
    parser.add_argument('outFile', type=str,  help='the output root file')
    args = parser.parse_args()
    
    readScan(args.inputFile,args.outFile)
