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

first_iter=1


"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/'
x_var='zero_thr'
std_index=9-1
n_iters=19
maximize='phi2'  # phi1, ph2, both
loc_ratios='lower right'  # posizione legenda plot rapporti
loc_deltas='lower right'  #  "          "      "   differenze
"""

"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/'
x_var='moma1_thr'
std_index=9-1
n_iters=20
maximize='phi2'  # phi1, ph2, both
loc_ratios='upper right'  # posizione legenda plot rapporti
loc_deltas='upper right'  #  "          "      "   differenze
"""

"""

base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/'
x_var='moma2_thr'
std_index=9-1
n_iters=20
maximize='phi2'  # phi1, ph2, both
loc_ratios='upper right'  # posizione legenda plot rapporti
loc_deltas='upper right'  #  "          "      "   differenze
"""



base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin/'
x_var='dmin'
std_index=8-1
n_iters=17
maximize='phi2'  # phi1, ph2, both
loc_ratios='lower left'  # posizione legenda plot rapporti
loc_deltas='lower left'  #  "          "      "   differenze


"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmax/'
x_var='dmax'
std_index=10-1-1
n_iters=21
maximize='phi2'  # phi1, ph2, both
loc_ratios='lower right'  # posizione legenda plot rapporti
loc_deltas='lower right'  #  "          "      "   differenze
"""


"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/'
x_var='weight_scale'
std_index=5-1-1 # il scondo -1 perche' escludo il primo punto
n_iters=20
maximize='phi2'  # phi1, ph2, both
loc_ratios='lower right'  # posizione legenda plot rapporti
loc_deltas='lower right'  #  "          "      "   differenze
first_iter=2  # escludo il primo punto!!
"""



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




def plot_all(baseRec1,outdir,folder):

    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    
    ##################
    # create plots
    ##################
    
    
    x=baseRec1.dict_rec[x_var]
        
    best_thr=x[baseRec1.best_index[0]]
    min_thr=x[baseRec1.min_index]
    max_thr=x[baseRec1.max_index]
    std_value=x[std_index]
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
    figAll_phi2=plt.figure(1,figsize=(10,7) )
    figAll_phi2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all=plt.subplot(111)
    ax_all.set_title("mod. factor phi2")
    plt.xlabel(x_var)
    plt.errorbar(x,y,yerr=y_err, fmt='o--',label='E='+str(dict_energy[folder])+' KeV - phi2 ' )
    #plt.legend( bbox_to_anchor=(1.1, 1.11) )

    #ax_all2=plt.subplot(122)
    #ax_all2.set_title("mod. factor phi1")

    plt.errorbar(x,y2,yerr=y2_err, fmt='*--',label='E='+str(dict_energy[folder])+' KeV - phi1' )


    #plt.xlabel(x_var)
    plt.legend(bbox_to_anchor=(1.1, 1.11))


    outfilePlot=base_dir+'summary3_short2.png'
    print ("outFile png =",outfilePlot)
    figAll_phi2.savefig(outfilePlot)



    mod1_np=np.array(y)
    mod1Err_np=np.array(y_err)
    nom_mod=mod1_np[std_index]
    nom_mod_err= mod1Err_np[std_index]  

    mod2_np=np.array(y2)
    mod2Err_np=np.array(y2_err)
    nom_mod2=mod2_np[std_index]
    nom_mod2_err= mod2Err_np[std_index]
    energyLabel='E='+str(dict_energy[folder])+' KeV'

    # plot cumulativo (all E), con mod relativa ai valori nominali  

    figAll_rel=plt.figure(2,figsize=(10,7) )
    figAll_rel.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all=plt.subplot(111)
    ax_all.set_title("mod. factor phi2 / mod factor phi2 [standard_parameters] ")
    plt.xlabel(x_var)

    ratio=mod1_np/nom_mod
    ratioErr=(  (1./nom_mod**2)*(mod1Err_np**2)+((mod1_np/(nom_mod**2))**2)*(nom_mod_err**2) )**0.5
    #plt.errorbar(x,ratio,yerr=ratioErr, fmt='o--',label='E='+folder+' KeV' )
    plt.errorbar(x,ratio, fmt='o--', label=energyLabel)
    plt.legend(loc=loc_ratios)
    outfilePlot=base_dir+'summary_ratios.png'
    print ("outFile png =",outfilePlot)
    figAll_rel.savefig(outfilePlot)

     
 # plot cumulativo (all E), con mod relativa ai valori nominali  PHI1

    figAll_rel2=plt.figure(3,figsize=(10,7) )
    figAll_rel2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all2=plt.subplot(111)
    ax_all2.set_title("mod. factor phi1 / mod factor phi1 [standard_parameters] ")
    plt.xlabel(x_var)

    ratio2=mod2_np/nom_mod2
    #ratioErr=(  (1./nom_mod**2)*(mod1Err_np**2)+((mod1_np/(nom_mod**2))**2)*(nom_mod_err**2) )**0.5
    #plt.errorbar(x,ratio,yerr=ratioErr, fmt='o--',label='E='+folder+' KeV' )
    plt.errorbar(x,ratio2, fmt='o--',label=energyLabel)
    plt.legend(loc=loc_ratios)
   
    outfilePlot=base_dir+'summary_ratios_phi1.png'
    print ("outFile png =",outfilePlot)
    figAll_rel2.savefig(outfilePlot)

    ###########################
    # plot tutte le energie con delta mod - phi2

    figAll_delta=plt.figure(4,figsize=(10,7) )
    figAll_delta.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all3=plt.subplot(111)
    ax_all3.set_title("mod. factor phi2 - mod factor phi2 [standard_parameters] ")
    plt.xlabel(x_var)

    delta=mod1_np-nom_mod
    plt.errorbar(x,delta, fmt='o--',label= energyLabel)
    plt.legend(loc=loc_deltas)
   
    outfilePlot=base_dir+'summary_deltas.png'
    print ("outFile png =",outfilePlot)
    figAll_delta.savefig(outfilePlot)


    ######################3
    # plot tutte le energie con delta mod - phi1

    figAll_delta2=plt.figure(5,figsize=(10,7) )
    figAll_delta2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all4=plt.subplot(111)
    ax_all4.set_title("mod. factor phi1 - mod factor phi1 [standard_parameters] - phi1 ")
    plt.xlabel(x_var)
    delta2=mod2_np-nom_mod2
    plt.errorbar(x,delta2, fmt='o--',label=energyLabel )
    plt.legend(loc=loc_deltas)
    outfilePlot=base_dir+'summary_deltas_phi1.png'
    print ("outFile png =",outfilePlot)
    figAll_delta2.savefig(outfilePlot)

    
    ###########################################
    # plot singole energie
    ###########################################
    # update 29.4.2020 -> accorpo plots modulazione, aggiungo plot fase, tolgo plot rapporti risoluzine e n.eventi

    
    fig01=plt.figure(figsize=(10,7))
    fig01.suptitle(title, fontsize=14)
    fig01.subplots_adjust(left=0.09, right=0.9, top=0.91, bottom=0.06,hspace=0.350)
    # PLOT  mod_2 vs zero threshold
    ax01=plt.subplot(311)
    ax01.set_title('modulation')
    plt.errorbar(x,y,yerr=y_err, fmt='bo--',label="modulation phi2")
    plt.xlabel(x_var,x=0.9)
    plt.ylabel('modulation')

    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )

    ax01.add_patch(rect)
    plt.errorbar(x,y2,yerr=y2_err, fmt='ro--',label="modulation phi1") # modulazione phi1

    ax01.set_ylim(   (min(y)-y_err[0]) - ((min(y)-y_err[0]))* 0.02, (max(y)+y_err[0])+ (max(y)+y_err[0])*0.02  )

    
    plt.legend(ncol=2)
  
    #plt.xlabel(x_var)
    #plt.ylabel('modulation_phi1')

    # plot phase...
    ax02=plt.subplot(312)
    ax02.set_title('phase')
    plt.errorbar(x, baseRec1.dict_rec['phase1']  ,yerr=baseRec1.dict_rec['phase1_err'] , fmt='ro--',label="phase phi1") # modulazione phi1
    plt.errorbar(x, baseRec1.dict_rec['phase2']  ,yerr=baseRec1.dict_rec['phase2_err'] , fmt='bo--',label="phase phi2") # modulazione phi1
  
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    plt.legend(ncol=2)
    plt.xlabel(x_var,x=0.9)
    plt.ylabel('phase [deg]')
    

    

    # plot chi2
    ax03=plt.subplot(313)
    ax03.set_title(r'modulation $\chi^2$')
    plt.plot(x, baseRec1.dict_rec['chi2_1'], 'ro--',label="chi2 phi1") # modulazione phi1
    plt.plot(x, baseRec1.dict_rec['chi2_2'], 'bo--',label="chi2 phi2") # modulazione phi1
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    
    plt.legend(ncol=2)
    plt.xlabel(x_var,x=0.9)
    plt.ylabel('Chi2')

    """
    ax04=plt.subplot(234)
    ax04.set_title('resolution')
    plt.errorbar(x, baseRec1.dict_rec['resolution2'], yerr= baseRec1.dict_rec['resolution2_err'], fmt='ro--',label="chi2 phi1")
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    
    plt.xlabel(x_var)
    plt.ylabel('PHA resolution FWHM')
    plt.legend()
 
    ax05=plt.subplot(235)
    ax05.set_title('ratio n_physical/n_raw')
    y=np.array(baseRec1.dict_rec['n_physical'])/np.array(baseRec1.dict_rec['n_raw'])
    plt.xlabel(x_var)
    plt.ylabel('fraction')
    plt.errorbar(x, y,fmt='ro--',label="n_physical/n_raw")
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    plt.legend()
    
    ax06=plt.subplot(236)
    ax06.set_title('ratio n_ecut/n_physical')
    y=np.array(baseRec1.dict_rec['n_ecut'])/np.array(baseRec1.dict_rec['n_physical'])
    plt.errorbar(x, y,fmt='ro--',label="n_ecut/n_physical")
    
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    plt.xlabel(x_var)
    plt.ylabel('fraction')
    plt.legend()
    """  


    
    outfilePlot=out_dir+'scan_summary_new.png'
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
 




#inizio scan su zero_thr
for folder in dirs:

    out_dir=base_dir+folder+'/'
    baseRec1=base_rec()
    
    for  i in range (first_iter,n_iters+1): #!!!!!!!!!!!!!!!!!!!!!!1
            
        work_dir=out_dir+str(i)
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
    
    print ('best x= ',baseRec1.dict_rec[x_var][baseRec1.best_index[0]], '  best mod2= ',baseRec1.dict_rec['modulation2'][baseRec1.best_index[0]])

    print ('std x=', baseRec1.dict_rec[x_var][std_index],  '  std mod2= ',baseRec1.dict_rec['modulation2'][std_index])

    print('\n x = ',baseRec1.dict_rec[x_var])
    print('y = ', baseRec1.dict_rec['modulation2'])


    
    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    plot_all(baseRec1,out_dir,folder) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    #fill variables for final plots (vs energy)
    energy.append(dict_energy[folder])
    best_val.append (baseRec1.dict_rec[x_var][baseRec1.best_index[0]] )
    best_val_up.append (baseRec1.dict_rec[x_var][baseRec1.max_index] )
    best_val_low.append (baseRec1.dict_rec[x_var][baseRec1.min_index] )

    
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][baseRec1.best_index[0]])
    n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][baseRec1.best_index[0]])
    mod1.append(baseRec1.dict_rec['modulation1'][baseRec1.best_index[0]])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][baseRec1.best_index[0]])
    mod2.append(baseRec1.dict_rec['modulation2'][baseRec1.best_index[0]])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][baseRec1.best_index[0]])
    phase1.append(baseRec1.dict_rec['phase1'][baseRec1.best_index[0]])
    phase1_err.append(baseRec1.dict_rec['phase1_err'][baseRec1.best_index[0]])
    phase2.append(baseRec1.dict_rec['phase2'][baseRec1.best_index[0]])
    phase2_err.append(baseRec1.dict_rec['phase2_err'][baseRec1.best_index[0]])
    


    resolution_opt.append(baseRec1.dict_rec['resolution2'][baseRec1.best_index[0]])
    n_final_opt.append(baseRec1.dict_rec['n_final'][baseRec1.best_index[0]])
    n_phys_opt.append(baseRec1.dict_rec['n_physical'][baseRec1.best_index[0]])

    

    
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
    std_val=baseRec1.dict_rec[x_var][std_index]



   

    
print("e=",energy)
print("best_thr=",best_val)
print('test string=',' '.join(map(str,energy)) )

###############
# write outFile.txt

outFileName=base_dir+'outScan_'+x_var+'_v2.txt'
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


outFile.write( ' '.join(map(str,best_val))+'\n'  )
outFile.write( ' '.join(map(str,best_val_up))+'\n'  )
outFile.write( ' '.join(map(str,best_val_low))+'\n'  )

outFile.close()




# final plots:

##########
#  FIG 1
##########

fig1=plt.figure(figsize=(20,10))
fig1.suptitle(x_var+" scan", fontsize=16)
fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
# n raw
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


    
outfilePlot=base_dir+'summary1_short2.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

#############
#  FIG 2
##############

fig2=plt.figure(figsize=(20,10))
fig2.suptitle(x_var+" scan", fontsize=16)
fig2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

# modulazione
ax1=plt.subplot(231)
ax1.set_title('modulation')
plt.errorbar(energy,mod2,yerr=mod2_err, fmt='ro--',label='phi 2 opt')
plt.errorbar(energy,mod2std,yerr=mod2std_err, fmt='bo--',label='phi 2 standard')
plt.errorbar(energy,mod1,yerr=mod1_err, fmt='co--',label='phi 1 opt ')
plt.errorbar(energy,mod1std,yerr=mod1std_err, fmt='go--',label='phi 1 standard ' )

plt.xlabel('energy [KeV]')
plt.ylabel('modulation')
plt.legend(loc='upper left')


#ratio modulation
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

ratio1=mod1_np/mod1std_np
ratio1Err=(  (1./mod1std_np**2)*(mod1Err_np**2)+((mod1_np/(mod1std_np**2))**2)*(mod1stdErr_np**2) )**0.5

ratio2=mod2_np/mod2std_np
ratio2Err=(  (1./mod2std_np**2)*(mod2Err_np**2)+((mod2_np/(mod2std_np**2))**2)*(mod2stdErr_np**2) )**0.5



ax2=plt.subplot(234)
#plt.errorbar(energy, np.array(mod1)/np.array(mod1std),  fmt='bo--',label='phi1  ratio opt/std ')
#plt.errorbar(energy, np.array(mod2)/np.array(mod2std),  fmt='ro--',label='phi2  ratio opt/std ')
plt.errorbar(energy, ratio1, yerr=ratio1Err,  fmt='bo--',label='phi1  ratio opt/std ')
plt.errorbar(energy, ratio2, yerr=ratio2Err,  fmt='ro--',label='phi2  ratio opt/std ')
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


ax4=plt.subplot(233)
ax4.set_title('best '+x_var)
plt.errorbar(energy,best_val, fmt='bo',label='best value')
ax4.fill_between(energy,best_val_low, best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax4.axhline(y=std_val,label='standard_value', linestyle='--',alpha=0.5)
plt.xlabel('energy [KeV]')
plt.ylabel('best '+x_var)
plt.legend()



# plot Delta phase

ax5=plt.subplot(236) #pippo
ax5.set_title('phase difference '+x_var)

deltaPhase1=phase1_np-phase1std_np
deltaPhase2=phase2_np-phase2std_np

plt.errorbar(energy,deltaPhase1, fmt='ro',label='delta phase 1')
plt.errorbar(energy,deltaPhase2, fmt='bo',label='delta phase 2')


plt.xlabel('energy [KeV]')
plt.ylabel('phase_opt - phase_std')
plt.legend()






outfilePlot=base_dir+'summary2_short2.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

                 
plt.show()
