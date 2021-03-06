
import ROOT
import numpy
from math import *
from array import array


def media_mobile (h1):        
    #nbins=h1.GetNbinsX()
    hF=h1.Clone()
    h2=h1.Clone()
    #h2.Rebin(2)
    
    hF.Reset()
    
    nbins=h2.GetNbinsX()
    
    for i in range (2,nbins-2):
         #ave=  (h1.GetBinContent(i)+ h1.GetBinContent(i-1)+h1.GetBinContent(i+1))/3.
         # ave=  (h1.GetBinContent(i)+  h1.GetBinContent(i-2)+  h1.GetBinContent(i-1)+h1.GetBinContent(i+1)+h1.GetBinContent(i+2)    )/5.

          ave=0
          for jj  in range (-1,1):
              ave=ave+h2.GetBinContent(i+jj)

          ave=ave/3.    
              
          hF.SetBinContent(i,ave)

    moda=hF.GetBinCenter(hF.GetMaximumBin())     
    gauss2=ROOT.TF1("gaus2","gaus",0,1000)
    gauss2.SetParameter(1,moda)
    gauss2.SetParameter(2,25)

    gauss2.SetParLimits(1,moda-25,moda+25)

    hF.Fit('gaus2',"M","",moda-25,moda+35)


    gaussAve=gauss2.GetParameter(1)

    landau=ROOT.TF1("landau","landau",0,150)
    landau.SetParameter(1, moda)
    hF.Fit('landau',"MER")

    landauAve=landau.GetParameter(1)
    print "gaussAve =",gaussAve," landauAve= ",landauAve
    
    valuesAve=[gaussAve,landauAve]

    #hF.Draw()
    #valore = raw_input('continue?')
    return valuesAve      





def fit2Gaussians(h1):
    
    gauss=ROOT.TF1("gaus","gaus",0,1000)

    moda=h1.GetBinCenter(h1.GetMaximumBin())
    gauss.SetParameter(1,moda)
    gauss.SetParameter(2,h1.GetRMS())

    gauss.SetParLimits(1,moda-50,moda+50)
    
    h1.Fit('gaus',"M","",0,150)
    
    meanG=gauss.GetParameter(1)
    RMSG=gauss.GetParameter(2)
    #RMSG=20
   
    h1.Fit('gaus',"M","",meanG-2.*RMSG,meanG+1.*RMSG)
    meanG=gauss.GetParameter(1)

   
    return meanG




def doAll(infile,outFile):


    print ("infile  =",infile)

       
    nCol=300
    nRows=352    

    #nCol=10
    #nRows=10    

    inRootFile=ROOT.TFile(infile,"open")

    x1=-8.99
    x2=8.01
    y1=-8.01
    y2=8.01

    nbinsX=int( (x2-x1)/0.025);
    nbinsY=int((y2-y1)/0.025)
    nbinsXhex=int( (x2-x1)/0.028);
    nbinsYhex=int((y2-y1)/0.028)
   
    h2hex=ROOT.TH2Poly("h2poly","",nbinsXhex,x1,x2,nbinsYhex,y1,y2)
    h2hex.Honeycomb(-8.99,-8.01,0.028,nbinsXhex ,nbinsYhex);

    h2=ROOT.TH2F("h2","",nbinsX,x1,x2,nbinsY,y1,y2)
    #h2pixel=ROOT.TH2F("h2pixel","",352,0,352,300,0,300)
    #h2pixelN=ROOT.TH2F("h2pixelN","",352,0,352,300,0,300)
  
    h2pixel=ROOT.TH2F("h2pixel","",300,0,300,352,0,352)
    h2pixelN=ROOT.TH2F("h2pixelN","",300,0,300,352,0,352)
  
    hz=ROOT.TH1F("hz","hz",500,0,500)
    hzRMS=ROOT.TH1F("hzRMS","hzRMS",500,0,500)

    hz=ROOT.TH1F("hz","hz",500,0,500)
    hzGauss=ROOT.TH1F("hzGauss","hzGauss",500,0,500)
    hzCic=ROOT.TH1F("hzCic","hzCic",500,0,500)
    hzCic3=ROOT.TH1F("hzCic3","hzCic3",500,0,500)
  
    print ("creating out root file: ",outFile)
    outRootFile=ROOT.TFile(outFile,"recreate")
   
    treeSimo=ROOT.TTree("tree","tree")
    mean= array('d', [0.])
    RMS= array('d', [0.])
    meanG=array('d', [0.])
    meanCic=array('d', [0.])
    meanCic3=array('d', [0.])
    gaussAve=array('d', [0.])
    landauAve=array('d', [0.])
    
    
    entries=array('d', [0.])
    pix_x=array('i', [0])
    pix_y=array('i', [0])
   

    treeSimo.Branch('mean', mean, 'mean/D') 
    treeSimo.Branch('RMS', RMS, 'RMS/D') 
    treeSimo.Branch('meanG', meanG, 'meanG/D') 
    treeSimo.Branch('meanCic', meanCic, 'meanCic/D') 
    treeSimo.Branch('meanCic3', meanCic3, 'meanCic3/D') 
    treeSimo.Branch('gaussAve', gaussAve, 'gaussAve/D') 
    treeSimo.Branch('landauAve', landauAve, 'landauAve/D') 
   
    treeSimo.Branch('entries', entries, 'entries/D')
    treeSimo.Branch('pix_x', pix_x, 'pix_x/I')
    treeSimo.Branch('pix_y', pix_y, 'pix_y/I')
 
     
    #for col in range (10,60):
    for col in range (1,nCol):
        print "col=", col

        for row in range (1,nRows):
        #for row in range (10,60):

            #print "col=", col," row= ",row
            pix_x[0]=col
            pix_y[0]=row
            
            hName='hist_'+str(col)+'_'+str(row)
            print "nome = ",hName
            hist=inRootFile.Get(hName)
            hist.Rebin(10)
            hist.Scale(1./float(hist.GetEntries()) )


            valuesAve=media_mobile (hist)
            gaussAve[0]=valuesAve[0]
            landauAve[0]=valuesAve[1]
            

            
            exp=ROOT.TF1("exp","expo",50,200)
            hist.Fit("exp","MER")
            meanCic3[0]=(exp.GetX(0.03)+exp.GetX(0.025))/2.
            meanCic[0]=exp.GetX(0.03)
            
            print "cic = ",meanCic[0], "cic3= ",meanCic3[0]

            mean[0]=  hist.GetMean()
            RMS[0]= hist.GetRMS()
            entries[0]= hist.GetEntries()

            meanG[0]= fit2Gaussians(hist)     
            #hist.Write() !!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            #wordlPix=ixpeGrid.pixelToWorld(row,col)
            #h2hex.Fill(wordlPix.x(),wordlPix.y(),mean)
            #h2.Fill(wordlPix.x(),wordlPix.y(),mean)
            h2pixel.Fill(col,row,mean[0])
            #h2pixel.Fill(row,col,mean[0])
            h2pixelN.Fill(col,row, entries[0])
            
            treeSimo.Fill()

            if entries[0]>0 and col>5 and col<nCol-5 and row>5 and row<nRows-5 :
            
                hz.Fill(mean[0])
                hzRMS.Fill( hist.GetRMS())
                hzGauss.Fill(meanG[0])
                hzCic.Fill(meanCic[0])
                hzCic3.Fill(meanCic3[0])
               
                
                

            del hist
            
    print ("writing histos... ")
    
    hz.Write()
    hzRMS.Write()
    h2pixel.Write()
    h2pixelN.Write()
    hzGauss.Write()
    hzCic.Write()
    hzCic3.Write()

    treeSimo.Write()
    print ("closing root file...")
    outRootFile.Close()
  
    return 0





if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('inputFile', type=str, help='the input fits file')
    parser.add_argument('outFile', type=str,  help='the output root file')
      
    args = parser.parse_args()

     
    doAll(args.inputFile,args.outFile)

    print ("THE END")

