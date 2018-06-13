
import ROOT
from array import array



class myTTree:    

    def __init__(self):

                
        self.event_id = array('i', [0])    
        self.num_tracks= array('i', [0])  # ancora da aggiungere!
        self.size = array('i', [0])
        self.pulse_height = array('i', [0])
        self.phi = array('d', [0.])
        self.absorptionX = array('d', [0.])
        self.absorptionY = array('d', [0.])
        self.barycenterX = array('d', [0.])
        self.barycenterY = array('d', [0.])
        self.mom2trans = array('d', [0.])
        self.mom2long = array('d', [0.])
        self.skweness = array('d', [0.])
        
        self.x_simo= array('d', [0.])
        self.y_simo=array('d', [0.])
        self.phiSimo=array('d', [0.])
        self.redChi2=array('d', [0.])
        self.sumpars=array('d',[0.]*7) # manca  il riempimento!!!!
    
        self.treeSimo=ROOT.TTree("t_recSimo","t_recSimo")

        
        
         
        self.treeSimo.Branch('event_id', self.event_id, 'event_id/I') 
        self.treeSimo.Branch('num_tracks', self.num_tracks, 'num_tracks/I')
        self.treeSimo.Branch('size', self.size, 'size/I')
        self.treeSimo.Branch('pulse_height',self.pulse_height, 'pulse_height/I')
        self.treeSimo.Branch('phi', self.phi, 'phi/D')
        self.treeSimo.Branch('absorptionX', self.absorptionX, 'absorptionX/D')
        self.treeSimo.Branch('absorptionY',self.absorptionY, 'absorptionY/D')
        self.treeSimo.Branch('barycenterY',self.barycenterY, 'barycenterY/D')
        self.treeSimo.Branch('barycenterX',self.barycenterX, 'barycenterX/D')
        self.treeSimo.Branch('mom2trans', self.mom2trans, 'mom2trans/D')
        self.treeSimo.Branch('mom2long', self.mom2long, 'mom2long/D')
        self.treeSimo.Branch('skweness', self.skweness, 'skweness/D')
        self.treeSimo.Branch('x_simo', self.x_simo, 'x_simo/D')
        self.treeSimo.Branch('y_simo', self.y_simo, 'y_simo/D')
        self.treeSimo.Branch('phiSimo', self.phiSimo, 'phiSimo/D')
        self.treeSimo.Branch('redChi2', self.redChi2, 'redChi2/D')
        self.treeSimo.Branch('sumpars', self.sumpars, 'sumpars[7]/D')

    def Fill(self,xpeSimo):

        self.event_id[0]=xpeSimo.event_id
        self.phi[0]=xpeSimo.phi1
        self.x_simo[0]=xpeSimo.xnew
        self.y_simo[0]=xpeSimo.ynew
        self.phiSimo[0]=xpeSimo.phiTang
        self.barycenterY[0]=xpeSimo.baricenter_Y
        self.barycenterX[0]=xpeSimo.baricenter_X

        self.absorptionX[0]=xpeSimo.conversion_point_X
        self.absorptionY[0]=xpeSimo.conversion_point_Y
        
        self.size[0]= xpeSimo.track.numHits()
        self.pulse_height[0]=xpeSimo.track.pulseHeight()
        #self.num_tracks=
               
        self.mom2trans[0]=xpeSimo.track.firstPassMomentsAnalysis().mom2trans()
        self.mom2long[0]=xpeSimo.track.firstPassMomentsAnalysis().mom2long()
        self.skweness[0]=xpeSimo.track.firstPassMomentsAnalysis().mom3long()
        self.redChi2[0]=xpeSimo.redChi2
        
        #self.sumpars=array('d',[0.]*7) # manca  il riempimento!!!!
        
        self.treeSimo.Fill()
        
