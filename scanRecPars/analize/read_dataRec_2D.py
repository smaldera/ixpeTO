from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl


from gpdswpy.binning import ixpeHistogram1d,ixpeHistogram2d

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

#import matplotlib.cm as cm
#from matplotlib.colors import Normalize



mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True


# PARAMETERS:

base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroMoma12/'
x_var='zero_thr'
var1='zero_thr'
var2='moma1_thr'
std_index=45-1
n_iters=138

maximize='both'  # phi1, ph2, both



dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV
#dirs=['001461','001436', '001416', '001388','001471','001361', '001333']
dirs=['001461','001436', '001416', '001361', '001333']


#dirs=['001461']



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
        self.min_x=1e20
        self.max_x=1e20
        self.min_y=1e20
        self.max_y=1e20


        

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







    def compute_indexes(self,maximum,mod,mod_err):          # mod e mod_err,sono numpy  array!

        x=np.array(baseRec1.dict_rec[var1])
        y=np.array(baseRec1.dict_rec[var2])
        
        index=np.where(mod==maximum)[0]
        maximum_err=mod_err[index[0]]     

        #masked_array=np.ma.masked_where( np.logical_and(z>maximum+maximum_err, z<maximum-maximum_err)  , z)
        masked_array=np.ma.masked_where( mod<maximum-maximum_err , mod)
    
        print ("mod=",mod)
        print("masked array = ",masked_array)
        mask=np.ma.getmask(masked_array)
        print("mask=",mask)

        masked_x=np.ma.array(x,mask=mask)
        masked_y=np.ma.array(y,mask=mask)
        
        print("masked_x=",masked_x)
        print ("x min=",np.min(masked_x)," max = ",np.max(masked_x))  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print("masked_y=",masked_y)
        print ("y min=",np.min(masked_y)," max = ",np.max(masked_y))  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              

        self.max_mod=maximum
        self.max_modErr=maximum_err
        self.best_index=index[0]
        self.min_x=np.min(masked_x)
        self.max_x= np.max(masked_x)
        self.min_y=np.min(masked_y)
        self.max_y= np.max(masked_y)
       

        
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
   # x_edges=np.array(range (14,30,2))
   # y_edges=np.array(range(14,39,2) )
    x_edges= np.linspace(4, 30,13)
    y_edges= np.linspace(4,40,18)


    x=np.array(baseRec1.dict_rec['zero_thr'])
    y=np.array(baseRec1.dict_rec['moma1_thr'])
    z=np.array(baseRec1.dict_rec['modulation2'])
    z_err=np.array(baseRec1.dict_rec['modulation2_err'])
    z1=np.array(baseRec1.dict_rec['modulation1'])
    z1_err=np.array(baseRec1.dict_rec['modulation1_err'])

    nraw=np.array(baseRec1.dict_rec['n_raw'])
    nphys=np.array(baseRec1.dict_rec['n_physical'])
    
    #mod1=np.array(z)
    rect=  patches.Rectangle(       (baseRec1.min_x, baseRec1.max_mod- baseRec1.max_modErr),   baseRec1.max_x-baseRec1.min_x,     2.*baseRec1.max_modErr,    facecolor='grey', edgecolor='none',alpha=0.3 )
    rect2=  patches.Rectangle(       (baseRec1.min_y, baseRec1.max_mod- baseRec1.max_modErr),   baseRec1.max_y-baseRec1.min_y,     2.*baseRec1.max_modErr,    facecolor='grey', edgecolor='none',alpha=0.3 )
     #                                 (x,y)                                                    width           height
    

    # creo tutti gli histos
    hist2d_mod2 = ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2,ztitle='mod phi_2')
    hist2d_mod1= ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2,ztitle='mod phi_1' )
    hist2d_chi2mod2= ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2, ztitle='chi2 phi_2' )
    hist2d_chi2mod1= ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2, ztitle='chi2 phi_1')
                   
    hist2d_res= ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2, ztitle='resolution')
    hist2d_nev= ixpeHistogram2d(x_edges, y_edges,  xtitle=var1, ytitle=var2,  ztitle='n_phys/n_raw' )
    # riempio histos
    hist2d_mod2.fill(x, y,weights=z)
    hist2d_mod1.fill(x, y,weights=z1)
    hist2d_chi2mod2.fill(x, y,weights=baseRec1.dict_rec['chi2_2'])
    hist2d_chi2mod1.fill(x, y,weights=baseRec1.dict_rec['chi2_1'])
    hist2d_res.fill(x, y,weights=baseRec1.dict_rec['resolution2'])
    hist2d_nev.fill(x, y,weights=nphys/nraw)

    
    
    #sample=(x,y)
    #data = np.vstack(sample).T
    #print ("DATA=",data.shape[0], "LEN wheights = ",len(z))

   
    #fig=plt.figure(figsize=(20,10))
    fig=plt.figure(figsize=(18,12))
    
    fig.suptitle(title, fontsize=16)
    fig.subplots_adjust(left=0.04,right=0.96, top=0.92, bottom=0.05,hspace=0.250,wspace=0.28)


    # PLOT mod phi2
    ax = fig.add_subplot(331)
   # hist2d_mod2.fill(x, y,weights=z)
    hist2d_mod2.plot(cmin=1e-10)
    #for i in range(len(y_edges)-1):
    #    for j in range(len(x_edges)-1):
    #        ax.text(x_edges[j]+0.5,y_edges[i]+0.5, hist2d_mod2[i,j], color="w", ha="center", va="center", fontweight="bold")

    
    if baseRec1.best_phi==2:
        ax.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='best values')
        ax.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values')
       
        rect3=  patches.Rectangle(       (baseRec1.min_x, baseRec1.min_y),   baseRec1.max_x-baseRec1.min_x,    baseRec1.max_y-baseRec1.min_y  ,    facecolor='grey', edgecolor='none',alpha=0.1 )
        #                                 (x,y)                                                    width                                height
        ax.add_patch(rect3)
    plt.legend(loc='lower right')
    
    # PLOT  mod_2 vs zero threshold
    ax2=plt.subplot(332)
    ax2.set_title('modulation phi_2')

    cmap = cm.autumn
    plt.errorbar(x,z,yerr=z_err,fmt='bo',label="modulation phi2")
    if baseRec1.best_phi==2:
        ax2.errorbar(x[baseRec1.best_index],z[baseRec1.best_index],yerr=z_err[baseRec1.best_index],marker='o', markersize=10, color='m', label='best_values')
        ax2.plot(x[std_index], z[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)      
        ax2.set_xlim(min(x)-1,max(x)+1)
        #plt.colorbar()
        ax2.add_patch(rect)

    plt.xlabel(var1)
    plt.ylabel('modulation_2')
    plt.legend(loc='lower right')
 
     # PLOT  mod_2 vs zero threshold
    ax3=plt.subplot(333)
    ax3.set_title('modulation phi_2')
    plt.errorbar(y,z,yerr=z_err, fmt='bo',label="modulation phi2")
    plt.xlabel(var2)
    plt.ylabel('modulation_2')
    if baseRec1.best_phi==2:
        ax3.errorbar(y[baseRec1.best_index],z[baseRec1.best_index],yerr=z_err[baseRec1.best_index],marker='o', markersize=10, color='m', label='max mod1')
        ax3.plot(y[std_index], z[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10 )
        ax3.set_xlim(min(y)-1,max(y)+1)   
        ax3.add_patch(rect2)
    plt.legend(loc='lower right')   
    #####
    # PLOT mod 1

    ax4 = fig.add_subplot(334)
    hist2d_mod1.plot(cmin=1e-10)
    if baseRec1.best_phi==1:
        
        ax4.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
        ax4.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
       
        rect3=  patches.Rectangle(       (baseRec1.min_x, baseRec1.min_y),   baseRec1.max_x-baseRec1.min_x,    baseRec1.max_y-baseRec1.min_y  ,    facecolor='grey', edgecolor='none',alpha=0.1 )
        #                                 (x,y)                                                    width                                height
        ax4.add_patch(rect3)
    plt.legend(loc='lower right')   
    # PLOT  mod_2 vs zero threshold
    ax5=plt.subplot(335)
 #   ax5.set_title('modulation phi_2')
    cmap = cm.autumn
    plt.errorbar(x,z1,yerr=z1_err,fmt='bo',label="modulation phi2")
    if baseRec1.best_phi==1:
        ax5.errorbar(x[baseRec1.best_index],z1[baseRec1.best_index],yerr=z1_err[baseRec1.best_index],marker='o', markersize=10, color='m', label='max mod1')
        ax5.plot(x[std_index], z1[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10 )
        ax5.set_xlim(min(x)-1,max(x)+1)
        #plt.colorbar()
        ax5.add_patch(rect)

    plt.xlabel(var1)
    plt.ylabel('modulation_2')
    plt.legend(loc='lower right')
     # PLOT  mod_2 vs zero threshold
    ax6=plt.subplot(336)
    #ax6.set_title('modulation phi_2')
    plt.errorbar(y,z1,yerr=z1_err, fmt='bo',label="modulation phi1")
    plt.xlabel(var2)
    plt.ylabel('modulation_2')
    if baseRec1.best_phi==1:
        ax6.errorbar(y[baseRec1.best_index],z1[baseRec1.best_index],yerr=z1_err[baseRec1.best_index],marker='o', markersize=10, color='m', label='max mod1')
        ax6.plot(y[std_index], z1[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
       
        ax6.set_xlim(min(y)-1,max(y)+1)
        ax6.add_patch(rect2)
    plt.legend(loc='lower right')    
   # riga 3
   # PLOT  chi2
    ax7=plt.subplot(337)
    plt.ylabel('modulation_2')
    if baseRec1.best_phi==1:
         hist2d_chi2mod1.plot(cmin=1e-10)
         
         #ax7.set_title('chi2  phi_1')
           
    if baseRec1.best_phi==2:
         hist2d_chi2mod2.plot(cmin=1e-10)
         #ax7.set_title('chi2  phi_2')
    ax7.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
    ax7.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
    plt.legend(loc='lower right')   
    # PLOT  resolution
    ax8=plt.subplot(338)    
    plt.ylabel('resolution FWHM')
    hist2d_res.plot(cmin=1e-10)
    ax8.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
    ax8.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
    plt.legend(loc='lower right')    
     
   # PLOT  n events
    ax9=plt.subplot(339)    
    plt.ylabel('n_raw/n_physical')
    hist2d_nev.plot(cmin=1e-10)
    ax9.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
    ax9.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
    plt.legend(loc='lower right')          
    
    outfilePlot=out_dir+'scan_summary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)
    
    



    
    
    

    ###
    # plot cumulativi...
    figAll_phi2=plt.figure(2,figsize=(20,10) )
    figAll_phi2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all=plt.subplot(221)
    ax_all.set_title("mod. factor phi2")
    plt.xlabel(var1)
    plt.errorbar(x,z,yerr=z_err, fmt='o',label='E='+str(dict_energy[folder])+' KeV' )
    plt.legend( bbox_to_anchor=(1.1, 1.11) )

    ax_all2=plt.subplot(222)
    ax_all2.set_title("mod. factor phi2")
    plt.errorbar(y,z,yerr=z_err, fmt='o',label='E='+str(dict_energy[folder])+' KeV' )
    plt.xlabel(var2)
    plt.legend(bbox_to_anchor=(1.1, 1.11))
   
    ax_all3=plt.subplot(223)
    ax_all3.set_title("mod. factor phi1")
    plt.xlabel(var1)
    plt.errorbar(x,z1,yerr=z1_err, fmt='o',label='E='+str(dict_energy[folder])+' KeV' )
    plt.legend( bbox_to_anchor=(1.1, 1.11) )

    ax_all4=plt.subplot(224)
    ax_all4.set_title("mod. factor phi1")
    plt.errorbar(y,z1,yerr=z1_err, fmt='o',label='E='+str(dict_energy[folder])+' KeV' )
    plt.xlabel(var2)
    plt.legend(bbox_to_anchor=(1.1, 1.11))
   # plt.show()


   

    
    outfilePlot=base_dir+'summary3.png'
    print ("outFile png =",outfilePlot)
    figAll_phi2.savefig(outfilePlot)

    """
    # plot singole energie
    
    fig01=plt.figure(figsize=(20,10))
    fig01.suptitle(title, fontsize=16)
    fig01.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    # PLOT  mod_2 vs zero threshold
    ax01=plt.subplot(231)
    ax01.set_title('modulation phi_2')
    plt.errorbar(x,y,yerr=y_err, fmt='bo--',label="modulation phi2")
    plt.xlabel(x_var)
    plt.ylabel('modulation_2')

    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
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
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    if baseRec1.best_phi==1:
        ax02.add_patch(rect)
    plt.legend(loc='lower right')
    plt.xlabel(x_var)
    plt.ylabel('modulation_phi1')

    ax03=plt.subplot(233)
    ax03.set_title(r'modulation $\chi^2$')
    plt.plot(x, baseRec1.dict_rec['chi2_1'], 'ro--',label="chi2 phi1") # modulazione phi1
    plt.plot(x, baseRec1.dict_rec['chi2_2'], 'bo--',label="chi2 phi2") # modulazione phi1
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    
    plt.legend()
    plt.xlabel(x_var)
    plt.ylabel('Chi2')


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
      


    
    outfilePlot=out_dir+'scan_summary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)
    #plt.show()

    """
   


    


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



best_var1=[]
best_var1_up=[]
best_var1_low=[]

best_var2=[]
best_var2_up=[]
best_var2_low=[]


n_raw_opt=[]
n_ecut_opt=[]

n_raw_std=[]
n_ecut_std=[]

resolution_opt=[]
resolution_std=[]




#inizio scan su zero_thr
for folder in dirs:

    out_dir=base_dir+folder+'/'
    baseRec1=base_rec()
    
    for  i in range (1,n_iters+1):
        work_dir=out_dir+str(i)
        file_out=work_dir+'/prova_out.txt'
        file_cfg=work_dir+'/config_simo.txt'
        
        print('reading file:',file_out)  
        baseRec1.read_file_rec(file_out,file_cfg)

  
    
      
    baseRec1.find_maxMod_index()
    print ('best_index=',baseRec1.best_index)
    print ('min_index=',baseRec1.min_index)
    print ('max_index=',baseRec1.max_index)
    print ('max_mod=',baseRec1.max_mod)
    print ('max_modErr=',baseRec1.max_modErr)
    print ('best_phi=',baseRec1.best_phi)


    print ('min_x=',baseRec1.min_x)
    print ('max_x=',baseRec1.max_x)
    print ('min_y=',baseRec1.min_y)
    print ('max_y=',baseRec1.max_y)
    
    
    
    print ('best x= ',baseRec1.dict_rec[x_var][baseRec1.best_index], '  best mod2= ',baseRec1.dict_rec['modulation2'][baseRec1.best_index])
    print ('std x=', baseRec1.dict_rec[x_var][std_index],  '  std mod2= ',baseRec1.dict_rec['modulation2'][std_index])


    

    
    #print('\n x = ',baseRec1.dict_rec[x_var])
    #print('y = ', baseRec1.dict_rec['modulation2'])


    
    title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    plot_all(baseRec1,out_dir,folder) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    #fill variables for final plots (vs energy)
    energy.append(dict_energy[folder])
    best_var1.append (baseRec1.dict_rec[var1][baseRec1.best_index] )
    best_var1_up.append (baseRec1.max_x )
    best_var1_low.append (baseRec1.min_x )

    best_var2.append (baseRec1.dict_rec[var2][baseRec1.best_index] )
    best_var2_up.append (baseRec1.max_y )
    best_var2_low.append (baseRec1.min_y )

    
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][baseRec1.best_index])
    n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][baseRec1.best_index])
    mod1.append(baseRec1.dict_rec['modulation1'][baseRec1.best_index])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][baseRec1.best_index])
    mod2.append(baseRec1.dict_rec['modulation2'][baseRec1.best_index])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][baseRec1.best_index])
    resolution_opt.append(baseRec1.dict_rec['resolution2'][baseRec1.best_index])

    

    mod1std.append(baseRec1.dict_rec['modulation1'][std_index])
    mod1std_err.append(baseRec1.dict_rec['modulation1_err'][std_index])
    mod2std.append(baseRec1.dict_rec['modulation2'][std_index])
    mod2std_err.append(baseRec1.dict_rec['modulation2_err'][std_index])
    resolution_std.append(baseRec1.dict_rec['resolution2'][std_index])
    n_raw_std.append(baseRec1.dict_rec['n_raw'][std_index])
    n_ecut_std.append(baseRec1.dict_rec['n_ecut'][std_index])

#print("e=",energy)
#print("best_thr=",best_zero_thr)





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
plt.errorbar(energy,best_var1, fmt='bo',label='best value')
ax11.fill_between(energy,best_var1_low, best_var1_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
plt.xlabel('energy [KeV]')
plt.ylabel('best '+var1)
plt.legend()
#plt.rc('grid',axes=True,  linestyle=":", color='grey')
#plt.grid(True,linestyle=':', color='grey')

# n raw
ax12=plt.subplot(212)
ax12.set_title('n events analized')
plt.errorbar(energy,n_raw_opt, fmt='bo--')
plt.xlabel('energy [KeV]')
plt.ylabel('n_raw')                

outfilePlot=base_dir+'summary1.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

#############
#  FIG 2
##############

fig2=plt.figure(figsize=(20,10))
fig2.suptitle(var1+" scan", fontsize=16)
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
ax4.set_title('Event number')
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

outfilePlot=base_dir+'summary2.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

                 
plt.show()


#plt.show()
 
