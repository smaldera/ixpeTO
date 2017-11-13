
import ROOT
import numpy
from math import *





def createHistogramsMatrix(nCol,nRows):
    pixHist= [[ROOT.TH1F]*nRows for i in range(nCol)]
    
    for col in range (0,nCol):
        for row in range (0,nRows):
            name='hist_'+str(col)+'_'+str(row)
            #print name
            pixHist[col][row]=ROOT.TH1F(name,name,500,0.,500.)
            ROOT.SetOwnership(pixHist[col][row], False) # importantissimo!!! se no il grabage collector python ci mette una vita!!!
    return pixHist


    



def doAll(infile,outFile):


    print ("infile  =",infile)

       
    nCol=300
    nRows=352    

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
    h2pixel=ROOT.TH2F("h2pixel","",300,0,300,352,0,352)

    hz=ROOT.TH1F("hz","hz",500,0,500)
    hzRMS=ROOT.TH1F("hzRMS","hzRMS",500,0,500)

    print ("creating out root file: ",outFile)
    outRootFile=ROOT.TFile(outFile,"recreate")
     
  
    
    for col in range (0,nCol):
        print "col=", col
        for row in range (0,nRows):

            #print "col=", col," row= ",row

            #hName=hist_299_342
            hName='hist_'+str(col)+'_'+str(row)
            #print "hname = ",hName
            hist=inRootFile.Get(hName)
            hist.Rebin(5)

            hist.Write()
            
            
            mean=  hist.GetMean()
            entries= hist.GetEntries()
            #wordlPix=ixpeGrid.pixelToWorld(row,col)
            #h2hex.Fill(wordlPix.x(),wordlPix.y(),mean)
            #h2.Fill(wordlPix.x(),wordlPix.y(),mean)
            h2pixel.Fill(col,row,mean)
        
            if entries>0 and col>5 and col<nCol-5 and row>5 and row<nRows-5 :
                hz.Fill(mean)
                hzRMS.Fill( hist.GetRMS())

            del hist
            
    print ("wring histos... ")
    
    hz.Write()
    hzRMS.Write()
    h2pixel.Write()
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

