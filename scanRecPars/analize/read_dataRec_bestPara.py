# loop su scan ixperecon 1dim (variando un solo parametro)
# produce i plot a energia fissa ( modulazione, risoluzione, chi2, etc al variare del parametro)
# produce i plot riassuntivi finali (parametro ottimale, modulazione best, etc) vs energia
# produce file riassuntivo di output

# add plot  fase


from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl

#import matplotlib.cm as cm
#from matplotlib.colors import Normalize



#mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True
#mpl.rcParams['font.size']=15  #!!!!!!!!!!!!!!!!!!!!!!!!!!


# PARAMETERS:

base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/'


#params_sets=['best3d','std']
#params_sets=['ixperecon_para','std']
params_sets=['ixpereconPara_LW','std']


index_best3d=0
std_index=1   # da ricavare dalla lista  params_sets (leggo i files out in successione)  !!!
outFileName='outBestParams_'+params_sets[0]+'_v2.txt' # aggiunto pha


dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV
dirs=['001461','001436', '001416', '001388','001471','001361', '001333']
#dirs=['001333']


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







    def compute_indexes(self,maximum,mod,mod_err):          # mod e mod_err sonu numpy  array!
        
        index=np.where(mod==maximum)[0]
        maximum_err=mod_err[index[0]]     
        index_interval=np.where( np.logical_and(mod<maximum+maximum_err, mod>maximum-maximum_err)  )    # DEVO USARE np.logical.and e non and  !!!!!
     
        
        min_index=np.amin(index_interval)
        max_index=np.amax(index_interval)

        self.max_mod=maximum
        self.max_modErr=maximum_err
        self.best_index=index
        self.min_index=min_index
        self.max_index=max_index

        
    def find_maxMod_index(self):

        mod2=np.array( self.modulation2)
        mod1=np.array( self.modulation1)

        max2=np.amax(mod2)
        max1=np.amax(mod1)

        if maximize=='phi2':
            mod2_err=np.array( self.modulation2_err)
            self.compute_indexes(max2,mod2,mod2_err)
            self.best_phi=2

        if maximize=='phi1':
                mod1_err=np.array( self.modulation1_err)
                self.compute_indexes(max1,mod1,mod1_err)
                self.best_phi=1
            
        
        if maximize=='both':
            if (max1>=max2):
                mod1_err=np.array( self.modulation1_err)
                self.compute_indexes(max1,mod1,mod1_err)
                self.best_phi=1

            else:
                mod2_err=np.array( self.modulation2_err)
                self.compute_indexes(max2,mod2,mod2_err)
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


phase1=[]
phase1_err=[]
phase2=[]
phase2_err=[]

phase1std=[]
phase1std_err=[]
phase2std=[]
phase2std_err=[]


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
 
pha_std=[]


############################################33
#inizio loop sulle folders
for folder in dirs:

    out_dir=base_dir+folder+'/'
    baseRec1=base_rec()


    #index=0
    for  par_set   in  params_sets: #!!!!!!!!!!!!!!!!!!!!!!1
            
        work_dir=out_dir+par_set
        file_out=work_dir+'/prova_out.txt'
        file_cfg=work_dir+'/config_simo.txt'
        
        print('reading file:',file_out)  
        baseRec1.read_file_rec(file_out,file_cfg)

        
        
     #   index+=1
        


    #baseRec1.find_maxMod_index()
       
    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 

    
    #fill variables for final plots (vs energy)
    energy.append(dict_energy[folder])

   # best_val.append (baseRec1.dict_rec[x_var][baseRec1.best_index[0]] )
   # best_val_up.append (baseRec1.dict_rec[x_var][baseRec1.max_index] )
   # best_val_low.append (baseRec1.dict_rec[x_var][baseRec1.min_index] )

    
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][index_best3d])
    n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][index_best3d])
    mod1.append(baseRec1.dict_rec['modulation1'][index_best3d])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][index_best3d])
    mod2.append(baseRec1.dict_rec['modulation2'][index_best3d])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][index_best3d])
    phase1.append(baseRec1.dict_rec['phase1'][index_best3d])
    phase1_err.append(baseRec1.dict_rec['phase1_err'][index_best3d])
    phase2.append(baseRec1.dict_rec['phase2'][index_best3d])
    phase2_err.append(baseRec1.dict_rec['phase2_err'][index_best3d])
    


    resolution_opt.append(baseRec1.dict_rec['resolution2'][index_best3d])
    n_final_opt.append(baseRec1.dict_rec['n_final'][index_best3d])
    n_phys_opt.append(baseRec1.dict_rec['n_physical'][index_best3d])

    

    
    mod1std.append(baseRec1.dict_rec['modulation1'][std_index])
    mod1std_err.append(baseRec1.dict_rec['modulation1_err'][std_index])
    mod2std.append(baseRec1.dict_rec['modulation2'][std_index])
    mod2std_err.append(baseRec1.dict_rec['modulation2_err'][std_index])

    phase1std.append(baseRec1.dict_rec['phase1'][std_index])
    phase1std_err.append(baseRec1.dict_rec['phase1_err'][std_index])
    phase2std.append(baseRec1.dict_rec['phase2'][std_index])
    phase2std_err.append(baseRec1.dict_rec['phase2_err'][std_index])
    
    resolution_std.append(baseRec1.dict_rec['resolution2'][std_index])
    n_raw_std.append(baseRec1.dict_rec['n_raw'][std_index])
    n_ecut_std.append(baseRec1.dict_rec['n_ecut'][std_index])
    n_final_std.append(baseRec1.dict_rec['n_final'][std_index])
    n_phys_std.append(baseRec1.dict_rec['n_physical'][std_index])

    pha_std.append(baseRec1.dict_rec['peak2'][std_index])



   

    
print("e=",energy)
print("best_thr=",best_val)
print('test string=',' '.join(map(str,energy)) )

###############
# write outFile.txt

#outFileName=base_dir+'outBestParams.txt'
outFileName=base_dir+outFileName #'outBestParams_v2.txt' # aggiunto pha
                   
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

outFile.write( ' '.join(map(str,phase1)) +'\n' )   # aggiunte informazioni phase!!
outFile.write( ' '.join(map(str,phase1_err))+'\n'  )
outFile.write( ' '.join(map(str,phase2))+'\n'  )
outFile.write( ' '.join(map(str,phase2_err))+'\n'  )
outFile.write( ' '.join(map(str,phase1std)) +'\n' )
outFile.write( ' '.join(map(str,phase1std_err))+'\n'  )
outFile.write( ' '.join(map(str,phase2std))+'\n'  )
outFile.write( ' '.join(map(str,phase2std_err))+'\n') 

outFile.write( ' '.join(map(str, pha_std))+'\n') 



                   
outFile.close()




# final plots:

##########
#  FIG 1
##########

fig1=plt.figure(figsize=(20,10))
fig1.suptitle("number of events", fontsize=16)
fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

ax11=plt.subplot(221)
ax11.set_title('n events analized')
plt.errorbar(energy,n_raw_opt, fmt='bo--')
plt.xlabel('energy [KeV]')
plt.ylabel('n_raw')                

ax12=plt.subplot(222)
ax12.set_title('n_ecut/n_raw')
plt.errorbar(energy,  np.array(n_ecut_opt)/np.array(n_raw_opt),  fmt='bo--',label='n_ecut/n_raw opt ')
plt.errorbar(energy, np.array(n_ecut_std)/np.array(n_raw_std),  fmt='ro--',label='n_ecut/n_raw  std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio n. events')
plt.legend()


ax13=plt.subplot(223)
ax13.set_title('n_final/n_physical')
plt.errorbar(energy,  ( np.array(n_final_opt)/np.array(n_phys_opt))    ,  fmt='bo--',label='ratio final/physical opt ')
plt.errorbar(energy,  ( np.array(n_final_std)/np.array(n_phys_std))    ,  fmt='ro--',label='ratio final/physical std')

plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()


ax14=plt.subplot(224)
plt.errorbar(energy,  ( np.array(n_ecut_opt)/np.array(n_raw_opt))/(np.array(n_ecut_std)/np.array(n_raw_std))    ,  fmt='bo--',label='event fraction  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()


    
outfilePlot=base_dir+'summary1.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

#############
#  FIG 2
##############

#creo numpy array dalle liste

mod1_np=np.array(mod1)
mod2_np=np.array(mod2)
mod1Err_np=np.array(mod1_err)
mod2Err_np=np.array(mod2_err)

mod1std_np=np.array(mod1std)
mod2std_np=np.array(mod2std)
mod1stdErr_np=np.array(mod1std_err)
mod2stdErr_np=np.array(mod2std_err)


phase1_np=np.array(phase1)
phase2_np=np.array(phase2)
phase1Err_np=np.array(phase1_err)
phase2Err_np=np.array(phase2_err)
phase1std_np=np.array(phase1std)
phase2std_np=np.array(phase2std)
phase1stdErr_np=np.array(phase1std_err)
phase2stdErr_np=np.array(phase2std_err)

ratio1=mod2_np/mod2std_np
ratio1Err=(  (1./mod1std_np**2)*(mod1Err_np**2)+((mod1_np/(mod1std_np**2))**2)*(mod1stdErr_np**2) )**0.5

#ratio2=mod2_np/mod2std_np
ratio2=mod2_np/mod1std_np

ratio2Err=(  (1./mod2std_np**2)*(mod2Err_np**2)+((mod2_np/(mod2std_np**2))**2)*(mod2stdErr_np**2) )**0.5


fig2=plt.figure(figsize=(20,10))
fig2.suptitle("mod ratios", fontsize=16)
fig2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

ax1=plt.subplot(111)
ax1.set_title('modulation')

#plt.errorbar(energy, np.array(mod1)/np.array(mod1std),  fmt='bo--',label='phi1  ratio opt/std ')
#plt.errorbar(energy, np.array(mod2)/np.array(mod2std),  fmt='ro--',label='phi2  ratio opt/std ')

#plt.errorbar(energy, ratio1, yerr=ratio1Err,  fmt='bo--',label='phi1  ratio best_3d/std ')
#plt.errorbar(energy, ratio2, yerr=ratio2Err,  fmt='ro--',label='phi2  ratio best_3d/std ')


plt.errorbar(energy, ratio1,  fmt='bo--',label='mod phi2 LW/phi2_std   ')
plt.errorbar(energy, ratio2,  fmt='ro--',label='mod  phi2_LW / phi1_std ')


plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()

outfilePlot=base_dir+'modRatios.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)


#plots resolution

# plot Delta phase

fig3=plt.figure(figsize=(10,7))
ax3=plt.subplot(111) #pippo
ax3.set_title('phase difference ')

deltaPhase1=phase1_np-phase1std_np
deltaPhase2=phase2_np-phase2std_np

plt.errorbar(energy,deltaPhase1, fmt='ro',label='delta phase 1')
plt.errorbar(energy,deltaPhase2, fmt='bo',label='delta phase 2')


plt.xlabel('energy [KeV]')
plt.ylabel('phase_opt - phase_std')
plt.legend()


outfilePlot=base_dir+'deletaPhase.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

                 
plt.show()
