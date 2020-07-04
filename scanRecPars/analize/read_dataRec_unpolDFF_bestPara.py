# loop su scan ixperecon 1dim (variando un solo parametro) - versione per MODULAZIONE SPURIA 
# produce i plot a energia fissa ( modulazione, risoluzione, chi2, etc al variare del parametro)
# produce i plot riassuntivi finali (parametro ottimale, modulazione best, etc) vs energia
# produce file riassuntivo di output



from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl


from readSummaryData import *


#import plotSummary.plotData1d



#import matplotlib.cm as cm
#from matplotlib.colors import Normalize



mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True
mpl.rcParams['font.size']=15  #!!!!!!!!!!!!!!!!!!!!!!!!!!



# PARAMETERS:

first_iter=1




base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/bestPara_phaParametrization/'
x_var='weight_scale'
std_index=0
#n_iters=20
maximize='phi2'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/outScan_weight_scale.txt'
loc_ratios='upper right'
loc_deltas='lower right'





energy_dirs=['2.04','2.29', '2.70', '2.98', '3.69', '5.89']
#energy_dirs=['5.89']
eps_dir={'2.04':['000647','000658']  ,'2.29':['000669','000677'], '2.70':['000686','000691'], '2.98':['000704', '000726'], '3.69':['000733','000740'], '5.89':['000744','000752']}   # dir angolo1 e angolo2









class base_rec():

    def __init__(self):
        
        self.peak2=[]
        self.peak2_err=[]
        self.resolution2=[]
        self.resolution2_err=[]
        self.phase1=[]
        self.phase1_err=[]
        self.modulation1=[]
        self.modulation1_err=[]
        self.chi2_1=[]
        self.phase2=[]
        self.phase2_err=[]
        self.modulation2=[]
        self.modulation2_err=[]
        self.chi2_2=[]
        self.n_raw=[]
        self.n_physical=[]
        self.n_ecut=[]
        self.n_final=[]
    
        self.zero_thr=[]
        self.moma1_thr=[]
        self.moma2_thr=[]
        self.dmin=[]
        self.dmax=[]
        self.weight_scale=[]

        
        self.dict_rec={'peak2':self.peak2,'peak2_err':self.peak2_err,'resolution2':self.resolution2,'resolution2_err':self.resolution2_err, 'phase1':self.phase1,'phase1_err':self.phase1_err,  'modulation1':self.modulation1, 'modulation1_err':self.modulation1_err,'chi2_1':self.chi2_1, 'phase2':self. phase2,'phase2_err':self.phase2_err, 'modulation2':self.modulation2, 'modulation2_err':self.modulation2_err,'chi2_2':self.chi2_2,'n_raw':self.n_raw, 'n_physical':self.n_physical,'n_ecut':self.n_ecut, 'n_final':self.n_final,'zero_thr':self.zero_thr,'moma1_thr':self.moma1_thr, 'moma2_thr':self.moma2_thr,'dmin':self.dmin,'dmax':self.dmax,'weight_scale':self.weight_scale}

        #variabili x analisi

        self.max_mod=-1
        self.max_modErr=-1
        self.best_index=-1
        self.min_index=-2
        self.max_index=-1
        self.best_phi=-1
        


        

    def read_file_rec(self,file_rec,file_cfg):

         f=open(file_rec,'r')
         vals=f.readlines()[0].split()
         keys=['peak2','peak2_err','resolution2','resolution2_err', 'phase1','phase1_err',  'modulation1', 'modulation1_err','chi2_1', 'phase2', 'phase2_err'  , 'modulation2', 'modulation2_err','chi2_2','n_raw', 'n_physical','n_ecut', 'n_final']  #!!!! specifico le keys a mano perche' per python <3.7 l'ordine non e' garantito!!!! (inoltre non riesco a iterare con un indice: dict.keys()[i] non gli piace?   )
        
         for i in range (0,len(keys)):
             self.dict_rec[keys[i]].append(float(vals[i]))

         f_cfg=open(file_cfg,'r')
         vals=f_cfg.readlines()[0].split()
         keys=['zero_thr','moma1_thr', 'moma2_thr','dmin','dmax','weight_scale']
                 
         for i in range (0,len(keys)):
             self.dict_rec[keys[i]].append(float(vals[2*i+1]))

             
    def cerca_contigui(self,index_list, best_index):
        """
        """ 
         
        i=np.where(index_list==best_index)[0][0]
     
        
        max_i=len(index_list)-1
        min_i=0
        if( i==len(index_list) or i==len(index_list)-1) :
             max_i=i
             print ("max_i messo a i")
        else:     
            for j in range (i,len(index_list)-1 ):
               delta=index_list[j+1]-index_list[j]
               if delta>1:
                    max_i=j
                    break

        if i==0 or i==1:
             min_i=0
        else:      
            for j1 in range (1,i):
                j2=i-j1
                delta=index_list[j2+1]-index_list[j2]
                if delta>1:
                    min_i=j2
                    break     

        print("cerca contigui: min_i=",min_i, "max_i = ",max_i)        
        return index_list[min_i:max_i+1]     # lo slice esclude l'estremo superiore!!!!     
      



    def compute_indexes(self,minimum,mod,mod_err):          # mod e mod_err sonu numpy  array!
    
        
        # QUA VOGLIO CERCARE IL MINIMO!!!!
        index=np.where(mod==minimum)[0]
        #maximum_err=mod_err[index[0]]
        min_err=mod_err[index[0]]
       
        index_interval=np.where( mod<minimum+min_err*(2**0.5)  )    # DEVO USARE np.logical.and e non and  !!!!!


        print ("type index_interval",type(index_interval))
        print ("index_interval=",index_interval[0])
        print ("len = ",len(index_interval[0]))

        index_intervalContiguo= self.cerca_contigui(index_interval[0], index)
        
        print ("index_intervalContiguo=",index_intervalContiguo)
        min_index=np.amin(index_intervalContiguo)
        max_index=np.amax(index_intervalContiguo)

        #min_index=np.amin(index_interval)
        #max_index=np.amax(index_interval)

        self.max_mod=minimum
        self.max_modErr=min_err
        self.best_index=index
        self.min_index=min_index
        self.max_index=max_index

        
    def find_maxMod_index(self):

        mod2=np.array( self.modulation2)
        mod1=np.array( self.modulation1)

        max2=np.amax(mod2)
        max1=np.amax(mod1)

        min1=np.amin(mod1)
        min2=np.amin(mod2)
        

        if maximize=='phi2':
            mod2_err=np.array( self.modulation2_err)
            #self.compute_indexes(max2,mod2,mod2_err)
            self.compute_indexes(min2,mod2,mod2_err)
            self.best_phi=2

        if maximize=='phi1':
                mod1_err=np.array( self.modulation1_err)
                self.compute_indexes(min1,mod1,mod1_err)
                self.best_phi=1
            
        
        if maximize=='both':
            if (min1<min2):
                mod1_err=np.array( self.modulation1_err)
                self.compute_indexes(min1,mod1,mod1_err)
                self.best_phi=1

            else:
                mod2_err=np.array( self.modulation2_err)
                self.compute_indexes(min2,mod2,mod2_err)
                self.best_phi=2
                  
              
        return 0   
        
        

##################################################





#liste per plot finali:
energy=[]

mod1=[]
mod1_err=[]
mod2=[]
mod2_err=[]

mod1std=[]
mod1std_err=[]
mod2std=[]
mod2std_err=[]



best_val=[]
best_val_up=[]
best_val_low=[]

n_raw_opt=[]
n_ecut_opt=[]

n_raw_std=[]
n_ecut_std=[]

resolution_opt=[]
resolution_std=[]

std_val=-100
n_final_std=[]
n_final_opt=[]
n_phys_std=[]
n_phys_opt=[]
 




#inizio scan 
for folder in energy_dirs:

    out_dir=base_dir+folder+'/'+eps_dir[folder][0]+'/'
    baseRec1=base_rec()
    
    #for  i in range (first_iter,n_iters+1): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    #for  i in range (2,n_iters+1):
       
    work_dir=out_dir
    file_out=work_dir+'/prova_out.txt'
    file_cfg=work_dir+'/config_simo.txt'
        
    print('reading file:',file_out)  
    baseRec1.read_file_rec(file_out,file_cfg)

    #print("dict_rec=",baseRec1.dict_rec)
      


    baseRec1.find_maxMod_index()
    print ('best_index=',baseRec1.best_index[0])
    print ('min_index=',baseRec1.min_index)
    print ('max_index=',baseRec1.max_index)
    print ('max_mod=',baseRec1.max_mod)
    print ('max_modErr=',baseRec1.max_modErr)
    print ('best_phi=',baseRec1.best_phi)
    

    
    title='energy= '+folder+' KeV' 


    #fill variables for final plots (vs energy)
    energy.append(float(folder))
    best_val.append (baseRec1.dict_rec[x_var][baseRec1.best_index[0]] )
    best_val_up.append (baseRec1.dict_rec[x_var][baseRec1.max_index] )
    best_val_low.append (baseRec1.dict_rec[x_var][baseRec1.min_index] )

    
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][baseRec1.best_index[0]])
    n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][baseRec1.best_index[0]])
    mod1.append(baseRec1.dict_rec['modulation1'][baseRec1.best_index[0]])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][baseRec1.best_index[0]])
    mod2.append(baseRec1.dict_rec['modulation2'][baseRec1.best_index[0]])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][baseRec1.best_index[0]])
    resolution_opt.append(baseRec1.dict_rec['resolution2'][baseRec1.best_index[0]])
    n_final_opt.append(baseRec1.dict_rec['n_final'][baseRec1.best_index[0]])
    n_phys_opt.append(baseRec1.dict_rec['n_physical'][baseRec1.best_index[0]])

    mod1std.append(baseRec1.dict_rec['modulation1'][std_index])
    mod1std_err.append(baseRec1.dict_rec['modulation1_err'][std_index])
    mod2std.append(baseRec1.dict_rec['modulation2'][std_index])
    mod2std_err.append(baseRec1.dict_rec['modulation2_err'][std_index])
    resolution_std.append(baseRec1.dict_rec['resolution2'][std_index])
    n_raw_std.append(baseRec1.dict_rec['n_raw'][std_index])
    n_ecut_std.append(baseRec1.dict_rec['n_ecut'][std_index])
    n_final_std.append(baseRec1.dict_rec['n_final'][std_index])
    n_phys_std.append(baseRec1.dict_rec['n_physical'][std_index])
    std_val=baseRec1.dict_rec[x_var][std_index]

    
   

    
print("e=",energy)
print("best_thr=",best_val)


###############
# write outFile.txt

outFileName=base_dir+'outBestPara_ixpereconParametrization.txt'
outFile=open(outFileName,'w')
outFile.write( ' '.join(map(str,energy))+'\n'  )
outFile.write( ' '.join(map(str,mod1)) +'\n' )
outFile.write( ' '.join(map(str,mod1_err))+'\n'  )
outFile.write( ' '.join(map(str,mod2))+'\n'  )
outFile.write( ' '.join(map(str,mod2_err))+'\n'  )
outFile.write( ' '.join(map(str,mod1std)) +'\n' )
outFile.write( ' '.join(map(str,mod1std_err))+'\n'  )
outFile.write( ' '.join(map(str,mod2std))+'\n'  )
outFile.write( ' '.join(map(str,mod2std_err))+'\n'  )

outFile.write( ' '.join(map(str,best_val))+'\n'  )
outFile.write( ' '.join(map(str,best_val_up))+'\n'  )
outFile.write( ' '.join(map(str,best_val_low))+'\n'  )

# modulazione spuria al paramtero ottimale per run polarizzati!
#outFile.write( ' '.join(map(str, mod_phi1_bestPol_list))+'\n'  )
#outFile.write( ' '.join(map(str, mod_phi2_bestPol_list))+'\n'  )
#outFile.write( ' '.join(map(str, mod_phi1_bestPolErr_list))+'\n'  )
#outFile.write( ' '.join(map(str, mod_phi2_bestPolErr_list))+'\n'  )
#outFile.write( ' '.join(map(str, bestParPol_list ))+'\n'  )




outFile.close()



