
import ROOT
import numpy

from gpdswswig.Io import ixpeLvl0bFitsFile
from gpdswswig.MonteCarlo import *
from gpdswswig.Recon import *
from gpdswswig.Event import ixpeEvent
from gpdswswig.Geometry import ixpeCartesianCoordinate
from math import *
from array import array




def doAll(inputFile, outFile, num_events):


   print "starting... infile=",inputFile," outFile = ",outFile," n_events= ",num_events        
           
   zeroSupThreshold=5
   soglia1=5
   soglia2=5
   clustering = ixpeClustering(zeroSupThreshold,6)

   #binary_file=ixpeLvl0bFitsFile("/home/maldera/IXPE/software/ixpeSw/gpdsw/bin/sim.fits")
   #binary_file=ixpeLvl0bFitsFile("/home/maldera/IXPE/data/MC/fits/sim_2.5KeV_-180deg.fits")
   binary_file=ixpeLvl0bFitsFile(inputFile)


   x_convMC = array('d', [0.])
   y_convMC= array('d', [0.])
   phiMC= array('d', [0.])
   EnergyMC= array('d', [0.])
   hpaMC= array('d', [0.])
   phiRec= array('d', [0.])
   xRec= array('d', [0.])
   yRec= array('d', [0.])
   dist2D= array('d', [0.])
   xRec_rt= array('d', [0.])
   yRec_rt= array('d', [0.])         
   raggioMin= array('d', [0.])         
   raggioMax= array('d', [0.])         

   tree_rec=ROOT.TTree("tree_rec","tree_rec")
   tree_rec.Branch('x_convMC',x_convMC , 'x_convMC/D')
   tree_rec.Branch('y_convMC',y_convMC , 'y_convMC/D')
   tree_rec.Branch('phiMC',phiMC , 'phiMC/D')
   tree_rec.Branch('EnergyMC',EnergyMC , 'EnergyMC/D')
   tree_rec.Branch('phiRec',phiRec , 'phiRec/D')
   tree_rec.Branch('xRec',xRec , 'xRec/D')
   tree_rec.Branch('yRec',yRec , 'yRec/D')
   tree_rec.Branch('xRec_rt',xRec_rt , 'xRec_rt/D')
   tree_rec.Branch('yRec_rt',yRec_rt , 'yRec_rt/D')
   tree_rec.Branch('dist2D', dist2D , 'dist2D/D')
   tree_rec.Branch('raggioMin', raggioMin , 'RaggioMin/D')
   tree_rec.Branch('raggioMax', raggioMax , 'RaggioMax/D')

   hDist=ROOT.TH1F("hDist","",1000,0,1)
   hDiffX=ROOT.TH1F("hDiffX","",2000,-1,1)
   hDiffY=ROOT.TH1F("hDiffY","",2000,-1,1)
   hPhi=ROOT.TH1F("hPhi","", 360, -3.14, 3.14)
   hPhiMC=ROOT.TH1F("hPhiMC","", 360,  -3.14,3.14  )
   hPhiDiff=ROOT.TH1F("hPhiDiff","", 360,  -3.14,3.14  )
   hDiffXrt=ROOT.TH1F("hDiffXrt","",2000,-1,1)
   hDiffYrt=ROOT.TH1F("hDiffYrt","",2000,-1,1)




   for i in xrange(num_events):
           try:    

               digiEvt = binary_file.next()    
           #except RuntimeError  as  e:
           except:
                #if str(e)=='StopIteration':
                break
                
           evt = ixpeEvent(digiEvt) 
           mcInfo= binary_file.readMcInfo(i+1)
           x_convMC[0]=mcInfo.absorbtionPointX
           y_convMC[0]=mcInfo.absorbtionPointY
           phiMC[0]=mcInfo.photoElectronPhi
           EnergyMC[0]=mcInfo.photonEnergy
           #hpa_mc[0]=1
           
           tracks = clustering.dbScan(evt)    
           track = tracks[0]
         
           #lancio Rec:
           #track.reconstruct(soglia1, soglia2, False)

           for kk in range(0,10):
                      raggioMin[0]=float(0.5*kk)

                      
                      for jj in range (0,14):
                                 raggioMax[0]=raggioMin[0]+float(0.5*jj)
                                 if raggioMax[0]> 8:
                                            break

                                 #print "raggioMin = ",raggioMin[0], "raggioMax = ",raggioMax[0]
                                 
                                 track.reconstruct(soglia1, soglia2, raggioMin[0], raggioMax[0],0.05, False)
 
                                 phiRec[0] =track.secondPassMomentsAnalysis().phi()
                                 xRec[0]= track.absorptionPoint().x()
                                 yRec[0]= track.absorptionPoint().y()

                                 dist2D[0]=sqrt( ((xRec[0]-x_convMC[0])**2)+((yRec[0]-y_convMC[0])**2))


                                 #calcolo x e y rec dopo aver rototraslato usndo x y mc e phiMC
      
                                 xRec_rt[0]= cos(phiMC[0])*(xRec[0]-x_convMC[0]) + sin(phiMC[0])*(yRec[0]-y_convMC[0])
                                 yRec_rt[0]= -sin(phiMC[0])*(xRec[0]-x_convMC[0]) + cos(phiMC[0])*(yRec[0]-y_convMC[0])

           
                                 hDist.Fill(dist2D[0])
                                 hDiffX.Fill(xRec[0]-x_convMC[0])
                                 hDiffY.Fill(yRec[0]-y_convMC[0])
                                 hPhi.Fill(phiRec[0])           
                                 hPhiMC.Fill(phiMC[0])
                                 hPhiDiff.Fill(phiRec[0]-phiMC[0])
                                 hDiffYrt.Fill(yRec_rt[0])
                                 hDiffXrt.Fill(xRec_rt[0])
         
                                 #print "McInfo.E =",mcInfo.photonEnergy, "x_convMC =",x_convMC," phiMC = ",phiMC," PhiRec = ",phiRec," x_rec =",xRec, "y_rec = ",yRec," dist2D = ",dist2D
                                 tree_rec.Fill()

           if i%1000 ==0:
                      print 'letti ',i,' eventi... '

   print 'creating output file: ',outFile                              
   outRootFile=ROOT.TFile(outFile,"recreate");

   hDist.Write()                      
   hDiffX.Write()                      
   hDiffY.Write()  
   hPhi.Write()
   hPhiMC.Write() 
   hPhiDiff.Write()
   hDiffXrt.Write()
   hDiffYrt.Write()


   tree_rec.Write()
   outRootFile.Close()
   print outFile, 'closed...'                  




if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('inputFile', type=str,  help='the input binary file (.fit)')
    parser.add_argument('outFile', type=str,  help='the input binary file (.root) [lo so che root non piace... ma io usa!!!] ')
    parser.add_argument('-n', '--num_events', type=int, default=100000, help = 'number of events to be processed')
    args = parser.parse_args()
    

    
    
    doAll(args.inputFile,args.outFile, args.num_events)
