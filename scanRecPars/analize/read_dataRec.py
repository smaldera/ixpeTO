from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl

#import matplotlib.cm as cm
#from matplotlib.colors import Normalize




mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True


#plt.rc('grid', linestyle=":", color='grey')



base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr/'
dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV
dirs=['001461','001436', '001416', '001388','001471','001361', '001333']







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

        zero_th=np.array( self.zero_thr)
        index=np.where(mod==maximum)[0]
        maximum_err=mod_err[index[0]]     
        index_interval=np.where( np.logical_and(mod<maximum+maximum_err, mod>maximum-maximum_err)  )    # DEVO USARE np.logical.and !!!!!
        min_index=np.amin(index_interval)
        max_index=np.amax(index_interval)

        self.max_mod=maximum
        self.max_modErr=maximum_err
        self.best_index=index
        self.min_index=min_index
        self.max_index=max_index

        
    def find_maxMod_index(self):

       # zero_th=np.array( self.zero_thr)
        mod2=np.array( self.modulation2)
        mod1=np.array( self.modulation1)

        max2=np.amax(mod2)
        max1=np.amax(mod1)


        if (max1>=max2):
            mod1_err=np.array( self.modulation1_err)
            self.compute_indexes(max1,mod1,mod1_err)
            self.best_phi=1

        else:
            mod2_err=np.array( self.modulation2_err)
            self.compute_indexes(max2,mod2,mod2_err)
            self.best_phi=2
            
            
                    
        #return index_final,  min_index,max_index, max_mod, max_nodErr,best_phi   
        return 0   
        
        

##################################################




def plot_all(baseRec1,outdir,folder):

    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    
    ##################
    # create plots
    ##################
    
    x=baseRec1.dict_rec['zero_thr']
    best_thr=x[baseRec1.best_index[0]]
    min_thr=x[baseRec1.min_index]
    max_thr=x[baseRec1.max_index]
    
    #print(x)
    
    print ('best_thr=',best_thr)
    print ('min_thr=',min_thr)
    print ('max_thr=',max_thr)
    
    
    y=baseRec1.dict_rec['modulation2']
    y_err=baseRec1.dict_rec['modulation2_err']

    y2=baseRec1.dict_rec['modulation1']
    y2_err=baseRec1.dict_rec['modulation1_err']
    
    
 
    rect=  patches.Rectangle(       (min_thr,baseRec1.max_mod- baseRec1.max_modErr),   max_thr-min_thr,     2.*baseRec1.max_modErr,    facecolor='grey', edgecolor='none',alpha=0.3 )
    #                                 (x,y)                                                    width           height

    ###
    # plot cumulativi...
    figAll_phi2=plt.figure(1,figsize=(20,10) )
    figAll_phi2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    #figAll_phi2.suptitle("mod. factor phi2 ", fontsize=16)
    ax_all=plt.subplot(121)
    ax_all.set_title("mod. factor phi2")
    plt.errorbar(x,y,yerr=y_err, fmt='o--',label='E='+str(dict_energy[folder])+' KeV' )
    plt.legend()
    ax_all2=plt.subplot(122)
    ax_all2.set_title("mod. factor phi1")
    plt.errorbar(x,y2,yerr=y2_err, fmt='o--',label='E='+str(dict_energy[folder])+' KeV' )
    plt.legend()

    
    
    fig01=plt.figure(figsize=(20,10))
    fig01.suptitle(title, fontsize=16)
    fig01.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    # PLOT  mod_2 vs zero threshold
    ax01=plt.subplot(231)
    ax01.set_title('modulation phi_2')
    plt.errorbar(x,y,yerr=y_err, fmt='bo--',label="modulation phi2")
    plt.xlabel('zero_suppression_th')
    plt.ylabel('modulation_2')

    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )

    if baseRec1.best_phi==2:
        ax01.add_patch(rect)
    plt.legend(loc='lower right')

    # PLOT  mod_1 vs zero threshold
    ax02=plt.subplot(232)
    ax02.set_title('modulation phi_1')
    plt.errorbar(x,y2,yerr=y2_err, fmt='ro--',label="modulation phi1") # modulazione phi1
    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    if baseRec1.best_phi==1:
        ax02.add_patch(rect)
    plt.legend(loc='lower right')
    plt.xlabel('zero_suppression_th')
    plt.ylabel('modulation_phi1')

    ax03=plt.subplot(233)
    ax03.set_title(r'modulation $\chi^2$')
    plt.plot(x, baseRec1.dict_rec['chi2_1'], 'ro--',label="chi2 phi1") # modulazione phi1
    plt.plot(x, baseRec1.dict_rec['chi2_2'], 'bo--',label="chi2 phi2") # modulazione phi1
    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    
    plt.legend()
    plt.xlabel('zero_suppression_th')
    plt.ylabel('Chi2')


    ax04=plt.subplot(234)
    ax04.set_title('resolution')
    plt.errorbar(x, baseRec1.dict_rec['resolution2'], yerr= baseRec1.dict_rec['resolution2_err'], fmt='ro--',label="chi2 phi1")
    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    
    plt.xlabel('zero_suppression_th')
    plt.ylabel('PHA resolution FWHM')
    plt.legend()
 
    ax05=plt.subplot(235)
    ax05.set_title('ratio n_physical/n_raw')
    y=np.array(baseRec1.dict_rec['n_physical'])/np.array(baseRec1.dict_rec['n_raw'])
    plt.xlabel('zero_suppression_th')
    plt.ylabel('fraction')
    plt.errorbar(x, y,fmt='ro--',label="n_physical/n_raw")
    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    plt.legend()
    
    ax06=plt.subplot(236)
    ax06.set_title('ratio n_ecut/n_physical')
    y=np.array(baseRec1.dict_rec['n_ecut'])/np.array(baseRec1.dict_rec['n_physical'])
    plt.errorbar(x, y,fmt='ro--',label="n_ecut/n_physical")
    
    plt.axvline(x=20,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    plt.xlabel('zero_suppression_th')
    plt.ylabel('fraction')
    plt.legend()
      


    
    outfilePlot=out_dir+'scan_summary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)
    #plt.show()





    


################################################


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



best_zero_thr=[]
best_zero_thr_up=[]
best_zero_thr_low=[]

n_raw_opt=[]
n_ecut_opt=[]

n_raw_std=[]
n_ecut_std=[]

resolution_opt=[]
resolution_std=[]


std_index=7

#inizio scan su zero_thr
for folder in dirs:

    out_dir=base_dir+folder+'/'
    baseRec1=base_rec()
    
    for  i in range (1,19):
        work_dir=out_dir+str(i)
        file_out=work_dir+'/prova_out.txt'
        file_cfg=work_dir+'/config_simo.txt'
        
        print('reading file:',file_out)  
        baseRec1.read_file_rec(file_out,file_cfg)

    #print("dict_rec=",baseRec1.dict_rec)
      


    baseRec1.find_maxMod_index()
    print ('best_index=',baseRec1.best_index)
    print ('min_index=',baseRec1.min_index)
    print ('max_index=',baseRec1.max_index)
    print ('max_mod=',baseRec1.max_mod)
    print ('max_modErr=',baseRec1.max_modErr)
    print ('best_phi=',baseRec1.best_phi)
        
    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    plot_all(baseRec1,out_dir,folder) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    #fill variables for final plots (vs energy)
    energy.append(dict_energy[folder])
    best_zero_thr.append (baseRec1.dict_rec['zero_thr'][baseRec1.best_index[0]] )
    best_zero_thr_up.append (baseRec1.dict_rec['zero_thr'][baseRec1.max_index] )
    best_zero_thr_low.append (baseRec1.dict_rec['zero_thr'][baseRec1.min_index] )
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][baseRec1.best_index[0]])
    n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][baseRec1.best_index[0]])
    mod1.append(baseRec1.dict_rec['modulation1'][baseRec1.best_index[0]])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][baseRec1.best_index[0]])
    mod2.append(baseRec1.dict_rec['modulation2'][baseRec1.best_index[0]])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][baseRec1.best_index[0]])
    resolution_opt.append(baseRec1.dict_rec['resolution2'][baseRec1.best_index[0]])

    

    mod1std.append(baseRec1.dict_rec['modulation1'][std_index])
    mod1std_err.append(baseRec1.dict_rec['modulation1_err'][std_index])
    mod2std.append(baseRec1.dict_rec['modulation2'][std_index])
    mod2std_err.append(baseRec1.dict_rec['modulation2_err'][std_index])
    resolution_std.append(baseRec1.dict_rec['resolution2'][std_index])
    n_raw_std.append(baseRec1.dict_rec['n_raw'][std_index])
    n_ecut_std.append(baseRec1.dict_rec['n_ecut'][std_index])

print("e=",energy)
print("best_thr=",best_zero_thr)





# final plots:

##########
#  FIG 1
##########

fig1=plt.figure(figsize=(20,10))
fig1.suptitle("zero_supp threshold scan", fontsize=16)
fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

# PLOT  mod_2 vs zero threshold
ax11=plt.subplot(211)
ax11.set_title('best zero supp threshod')
plt.errorbar(energy,best_zero_thr, fmt='bo',label='best value')
ax11.fill_between(energy,best_zero_thr_low, best_zero_thr_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
plt.xlabel('energy [KeV]')
plt.ylabel('best_threshold')
plt.legend()
#plt.rc('grid',axes=True,  linestyle=":", color='grey')
#plt.grid(True,linestyle=':', color='grey')

# n raw
ax12=plt.subplot(212)
ax12.set_title('n events analized')
plt.errorbar(energy,n_raw_opt, fmt='bo--')
plt.xlabel('energy [KeV]')
plt.ylabel('n_raw')                


#############
#  FIG 2
##############

fig2=plt.figure(figsize=(20,10))
fig2.suptitle("zero_supp threshold scan", fontsize=16)
fig2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

# modulazione
ax1=plt.subplot(231)
ax1.set_title('modulation')
plt.errorbar(energy,mod2,yerr=mod2_err, fmt='ro--',label='phi 2 opt')
plt.errorbar(energy,mod2std,yerr=mod2std_err, fmt='bo--',label='phi 2 standard')
plt.errorbar(energy,mod1,yerr=mod1_err, fmt='bo--',label='phi 1 opt ')
plt.errorbar(energy,mod1std,yerr=mod1std_err, fmt='go--',label='phi 1 standard ' )

plt.xlabel('energy [KeV]')
plt.ylabel('modulation')
plt.legend()

ax2=plt.subplot(234)
plt.errorbar(energy, np.array(mod1)/np.array(mod1std),  fmt='bo--',label='phi1  ratio opt/std ')
plt.errorbar(energy, np.array(mod2)/np.array(mod2std),  fmt='ro--',label='phi2  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()

#plots resolution

ax3=plt.subplot(232)
ax3.set_title('Resolution')
plt.errorbar(energy, resolution_opt,  fmt='bo--',label='resolution opt ')
plt.errorbar(energy, resolution_std,  fmt='ro--',label='resolution std ')
plt.xlabel('energy [KeV]')
plt.ylabel('resolution')
plt.legend()

plt.subplot(235)
plt.errorbar(energy, np.array(resolution_opt)/np.array(resolution_std),  fmt='bo--',label='resolution  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()

#plots n. of events

ax4=plt.subplot(233)
ax4.set_title('Resolution')
plt.errorbar(energy,  np.array(n_ecut_opt)/np.array(n_raw_opt),  fmt='bo--',label='n_ecut/n_raw opt ')
plt.errorbar(energy, np.array(n_ecut_std)/np.array(n_raw_std),  fmt='ro--',label='n_ecut/n_raw  std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio n. events')
plt.legend()

plt.subplot(236)
plt.errorbar(energy,  ( np.array(n_ecut_opt)/np.array(n_raw_opt))/(np.array(n_ecut_std)/np.array(n_raw_std))    ,  fmt='bo--',label='event fraction  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()



                 
plt.show()
