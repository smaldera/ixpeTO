
import ROOT
import numpy
from math import *
from array import array

from astropy.io import fits
from gpdswpy.binning import ixpePixelMap

def doAll(args):

    #xedges=numpy.arange(300)
    ##yedges=numpy.arange(350)
    #tedges=numpy.arange(1)
    
    #print xedges
    #print yedges
    #print tedges

    inROOTfile=ROOT.TFile(args.inputFile,"open")
    h2=inROOTfile.Get(args.histName)  
   
    
        
    pix_map=ixpePixelMap.create_empty(300,352)
    
    for x in range (0,300):
          for y in range (0,352):
              global_bin = h2.GetBin(x,y);
              value=h2.GetBinContent(global_bin) 
              pix_map.fill(x,y,value)
    
    pix_map.normalize()
    pix_map.plot()         
 
    #write as fits!!!
    pix_map.write_fits(args.outFile)

    


if __name__ == '__main__':
    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('inputFile', type=str, help='the input root file')
    parser.add_argument('outFile', type=str,  help='out fits file')
    parser.add_argument('--histName', type=str, default='h2pixel',  help='histogram name')
    
    args = parser.parse_args()

     
    doAll(args)

    print ("THE END")

