
import ROOT
import numpy
from math import *


from gpdswswig.Io import ixpeLvl0bFitsFile
from gpdswswig.MonteCarlo import *
from gpdswswig.Recon import *
from gpdswswig.Event import ixpeEvent
from gpdswswig.Geometry import *



def createHistogramsMatrix(nCol,nRows):
    pixHist= [[ROOT.TH1F]*nRows for i in range(nCol)]
    
    for col in range (0,nCol):
        for row in range (0,nRows):
            name='hist_'+str(col)+'_'+str(row)
            #print name
            pixHist[col][row]=ROOT.TH1F(name,name,500,0.,500.)
            ROOT.SetOwnership(pixHist[col][row], False) # importantissimo!!! se no il grabage collector python ci mette una vita!!!
    return pixHist


    



def doAll(infileList,outFile,nEvent,zeroSupThreshold):


    print ("infileList  =",infileList)

    ixpeGrid= ixpeGeometrySvc.xpolAsicGrid()

    
    nCol=300
    nRows=352    
    pixHist= createHistogramsMatrix(nCol,nRows)

   
    
    clustering = ixpeClustering(zeroSupThreshold,5)

    for filename in infileList:
           
        print ("processing ",filename)
        binary_file=ixpeLvl0bFitsFile(filename)  
       
        for i in range (0, nEvent):
            try:
                digiEvt = binary_file.next() 
            except:
                #if str(e)=='StopIteration':
                break
            if i%10000==0:
                print "event = ",i
            evt= ixpeEvent(digiEvt)
            tracks = clustering.dbScan(evt)    
            track = tracks[0]

            # lancio rec...
            track.reconstruct(zeroSupThreshold,zeroSupThreshold, False) 

            absX=track.absorptionPoint().x()
            absY=track.absorptionPoint().y()


            hit=track.hits()
            n_hits=track.numHits()
           
            x=numpy.array([0.]*n_hits)
            y=numpy.array([0.]*n_hits)
            adc=numpy.array([0.]*n_hits)
   
      
            for i in range (0,n_hits):
                x[i]=hit[i].x
                y[i]=hit[i].y
                adc[i]=hit[i].pulseHeight
                dist=sqrt(((x[i]-absX)**2)+((y[i]-absY)**2))
                if  dist<0.05:

                    pixel=ixpeGrid.worldToPixel(x[i],y[i])
                    
                    col=pixel.column()
                    row=pixel.row()
                    pixHist[col][row].Fill(adc[i])
            # end for hits
        # end for events
       
       
        binary_file.close()
       
        
    #end loop files   

    print ("dopo loop files... " )
    
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
        for row in range (0,nRows):

            #print "col=", col," row= ",row
            
            pixHist[col][row].Write()
            mean=  pixHist[col][row].GetMean()
            entries= pixHist[col][row].GetEntries()
            wordlPix=ixpeGrid.pixelToWorld(row,col)
            h2hex.Fill(wordlPix.x(),wordlPix.y(),mean)
            h2.Fill(wordlPix.x(),wordlPix.y(),mean)
            h2pixel.Fill(col,row,mean)
        
            if entries>0:
                hz.Fill(mean)
                hzRMS.Fill( pixHist[col][row].GetRMS())

                   
    print ("wring histos... ")
    h2.Write()
    h2hex.Write()
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
    parser.add_argument('inputFile', type=str, nargs='+',   help='the input fits file')
    parser.add_argument('outFile', type=str,  help='the output root file')
    parser.add_argument('-n', '--num_events', type=int, default=10000000, help = 'number of events to be processed')
    parser.add_argument('-z', '--zero_suppression', type=int, default=5,help = 'zero-suppression threshold')
    
    args = parser.parse_args()

     
    doAll(args.inputFile,args.outFile,args.num_events,args.zero_suppression)

    print ("ora ho veramente finito...")

