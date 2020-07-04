#loop su scan ixperecon 2dim (variando 2 parametri)
# produce i plot a energia fissa ( modulazione, risoluzione, chi2, etc al variare dei parametri)
# produce i plot riassuntivi finali (parametri ottimal, modulazione best, etc) vs energia
# produce file riassuntivo di output



from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm

from gpdswpy.binning import ixpeHistogram1d,ixpeHistogram2d

import subprocess


mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True


# PARAMETERS:




base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/LWbins'
var1='dmin'
var2='weight_scale'
std_index=87-1
n_iters=132
maximize='phi2'  # phi1, ph2, both
x_edges= np.linspace(0., 2.2,12)     #dmin
y_edges= np.linspace(0., 0.24,13) #weights
var1_padding=0.1
var2_padding=0.05
LW_bins=24+1




dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV
dirs=['001461','001436', '001416', '001388','001471','001361', '001333']
#dirs=[ '001333']



n_iter_folders=0



def run_command(cmd):
         
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)



class base_rec():

    def __init__(self):
        
        #self.peak2=[]
        #self.peak2_err=[]
        #self.resolution2=[]
        #self.resolution2_err=[]




        self.pha_mean=[]
        self.pha_rms=[]
        
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
        #self.n_physical=[]
        #self.n_ecut=[]
        self.n_final=[]

        self.LWbin_start=[]
        self.LWbin_stop=[]
        self.LWmean=[]
        self.LWrms=[]
        
        self.zero_thr=[]
        self.moma1_thr=[]
        self.moma2_thr=[]
        self.dmin=[]
        self.dmax=[]
        self.weight_scale=[]

        
        self.dict_rec={'pha_mean':self.pha_mean,'pha_rms':self.pha_rms, 'phase1':self.phase1,'phase1_err':self.phase1_err,  'modulation1':self.modulation1, 'modulation1_err':self.modulation1_err,'chi2_1':self.chi2_1, 'phase2':self. phase2,'phase2_err':self.phase2_err, 'modulation2':self.modulation2, 'modulation2_err':self.modulation2_err,'chi2_2':self.chi2_2,'n_raw':self.n_raw, 'n_final':self.n_final, 'LWbin_start':self.LWbin_start,'LWbin_stop':self.LWbin_stop,'LWmean':self.LWmean,'LWrms':self.LWrms,'zero_thr':self.zero_thr,'moma1_thr':self.moma1_thr, 'moma2_thr':self.moma2_thr,'dmin':self.dmin,'dmax':self.dmax,'weight_scale':self.weight_scale}

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


         #print ("=====> open file:",file_rec)
         f=open(file_rec,'r')
         vals=f.readlines()[0].split()
         keys=[ 'pha_mean','pha_rms', 'phase1','phase1_err',  'modulation1', 'modulation1_err','chi2_1', 'phase2', 'phase2_err'  , 'modulation2', 'modulation2_err','chi2_2', 'n_final',  'n_raw', 'LWbin_start','LWbin_stop','LWmean','LWrms']  #!!!! specifico le keys a mano perche' per python <3.7 l'ordine non e' garantito!!!! (inoltre non riesco a iterare con un indice: dict.keys()[i] non gli piace?   )
        
         for i in range (0,len(keys)):
             self.dict_rec[keys[i]].append(float(vals[i]))



         #print("====> open_file",file_cfg)    
         f_cfg=open(file_cfg,'r')
         vals=f_cfg.readlines()[0].split()
         keys=['zero_thr','moma1_thr', 'moma2_thr','dmin','dmax','weight_scale']
                 
         for i in range (0,len(keys)):
             self.dict_rec[keys[i]].append(float(vals[2*i+1])) # 2*1+1 perche' ci sono anche le descrizioni delle variabili che "salto"







    def compute_indexes(self,maximum,mod,mod_err):          # mod e mod_err,sono numpy  array!

        x=np.array(baseRec1.dict_rec[var1])
        y=np.array(baseRec1.dict_rec[var2])
        
        index=np.where(mod==maximum)[0]
        maximum_err=mod_err[index[0]]     

        #masked_array=np.ma.masked_where( np.logical_and(z>maximum+maximum_err, z<maximum-maximum_err)  , z)
        masked_array=np.ma.masked_where( mod<maximum-maximum_err , mod)
    
       # print ("mod=",mod)
       # print("masked array = ",masked_array)
        mask=np.ma.getmask(masked_array)
        #print("mask=",mask)

        masked_x=np.ma.array(x,mask=mask)
        masked_y=np.ma.array(y,mask=mask)
        
       # print("masked_x=",masked_x)
       # print ("x min=",np.min(masked_x)," max = ",np.max(masked_x))  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       # print("masked_y=",masked_y)
       # print ("y min=",np.min(masked_y)," max = ",np.max(masked_y))  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              

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




def plot_all(baseRec1,outdir,LWfolder):

    #title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    title=LWfolder 
  
    ##################
    # create plots
    ##################
   # x_edges=np.array(range (14,30,2))
   # y_edges=np.array(range(14,39,2) )
#    x_edges= np.linspace(4, 30,13)
#    y_edges= np.linspace(4,40,18)


    x=np.array(baseRec1.dict_rec[var1])
    y=np.array(baseRec1.dict_rec[var2])
    z=np.array(baseRec1.dict_rec['modulation2'])
    z_err=np.array(baseRec1.dict_rec['modulation2_err'])
    z1=np.array(baseRec1.dict_rec['modulation1'])
    z1_err=np.array(baseRec1.dict_rec['modulation1_err'])

    nraw=np.array(baseRec1.dict_rec['n_raw'])
    #nphys=np.array(baseRec1.dict_rec['n_physical'])
    
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

    #hist2d_res.fill(x, y,weights=baseRec1.dict_rec['resolution2'])
    #hist2d_nev.fill(x, y,weights=nphys/nraw)

    
    
    #sample=(x,y)
    #data = np.vstack(sample).T
    #print ("DATA=",data.shape[0], "LEN wheights = ",len(z))

   
    #fig=plt.figure(figsize=(20,10))
    fig=plt.figure(figsize=(18,12))
    
    fig.suptitle(title, fontsize=16)
    fig.subplots_adjust(left=0.04,right=0.96, top=0.92, bottom=0.05,hspace=0.250,wspace=0.28)


    # PLOT mod phi2
    ax = fig.add_subplot(331)
  
    hist2d_mod2.plot(cmin=1e-10)
    # tentativo di scrivere il  valore del bin... bisogna capire come prendere i valori dall'ixpe2dhisto
    #for i in range(len(y_edges)-1):
    #    for j in range(len(x_edges)-1):
    #        ax.text(x_edges[j]+0.5,y_edges[i]+0.5, hist2d_mod2[i,j], color="w", ha="center", va="center", fontweight="bold")

    
    if baseRec1.best_phi==2:
        ax.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='best values')
        ax.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values')
       
        rect3=  patches.Rectangle(       (baseRec1.min_x, baseRec1.min_y),   baseRec1.max_x-baseRec1.min_x,    baseRec1.max_y-baseRec1.min_y  ,    facecolor='grey', edgecolor='none',alpha=0.1 )
        #                                 (x,y)                                                    width                                height
        ax.add_patch(rect3)
    #plt.legend(loc='lower right')
    plt.legend(loc='upper left')
   
    
    # PLOT  mod_2 vs zero threshold
    ax2=plt.subplot(332)
    ax2.set_title('modulation phi_2')

    cmap = cm.autumn
    plt.errorbar(x,z,yerr=z_err,fmt='bo',label="modulation phi2")
    if baseRec1.best_phi==2:
        ax2.errorbar(x[baseRec1.best_index],z[baseRec1.best_index],yerr=z_err[baseRec1.best_index],marker='o', markersize=10, color='m', label='best_values')
        ax2.plot(x[std_index], z[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)      
        ax2.set_xlim(min(x)-var1_padding,max(x)+var2_padding)
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
        ax3.set_xlim(min(y)-var2_padding,max(y)+var2_padding)   
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
        ax5.set_xlim(min(x)-var1_padding,max(x)+var1_padding)
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
       
        ax6.set_xlim(min(y)-var2_padding,max(y)+var2_padding)
        ax6.add_patch(rect2)
    plt.legend(loc='lower right')    
   # riga 3
   # PLOT  chi2
    ax7=plt.subplot(337)
    plt.ylabel('modulation_2')
    if baseRec1.best_phi==1:
         hist2d_chi2mod1.plot(cmin=1e-10,cmax=10)
         
         #ax7.set_title('chi2  phi_1')
           
    if baseRec1.best_phi==2:
         hist2d_chi2mod2.plot(cmin=1e-10,cmax=10)
         #ax7.set_title('chi2  phi_2')
    ax7.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
    ax7.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
    plt.legend(loc='lower right')   

    ## PLOT  resolution
    #ax8=plt.subplot(338)    
    #plt.ylabel('resolution FWHM')
    #hist2d_res.plot(cmin=1e-10)
    #ax8.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
    #ax8.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
    #plt.legend(loc='lower right')    
     
   # PLOT  n events
#    ax9=plt.subplot(339)    
#    plt.ylabel('n_raw/n_physical')
#    hist2d_nev.plot(cmin=1e-10)
#    ax9.plot(x[baseRec1.best_index], y[baseRec1.best_index], marker='o', markersize=10, color='m', label='max mod1')
#    ax9.plot(x[std_index], y[std_index], marker='o', markersize=10, color='c', label='std values', zorder=10)
#    plt.legend(loc='lower right')          
    #outfilePlot=out_dir+'scan_summary.png'

    outfilePlot=outdir+'scan_summary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)
    
    



    
    
    
    std_value1=x[std_index]
    std_value2=y[std_index]

    
   
    ###
    # plot cumulativi...
    leg_x=1.27
    leg_y=1.1

    figAll_phi2=plt.figure(2,figsize=(20,10) )
    figAll_phi2.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.09,hspace=0.250,wspace=0.45)
    ax_all=plt.subplot(221)
    ax_all.set_title("mod. factor phi2 vs "+var1)
    plt.xlabel(var1)
    plt.errorbar(x,z,yerr=z_err, fmt='o',label=LW_folder )
    plt.plot(x[baseRec1.best_index], z[baseRec1.best_index], marker='o', markersize=10,mfc='none', color='k', zorder=10)
       
    if n_iter_folders==1:
        plt.axvline(x=std_value1,label='standard_value', linestyle='--',alpha=0.5)
        plt.plot(x[baseRec1.best_index], z[baseRec1.best_index], marker='o', markersize=10,mfc='none', color='k', label='best values',zorder=10)
        
    plt.legend( bbox_to_anchor=(leg_x, leg_y ) )
    
    
    ax_all2=plt.subplot(222)
    ax_all2.set_title("mod. factor phi2 vs "+var2)
    plt.errorbar(y,z,yerr=z_err, fmt='o',label=LW_folder )
    plt.xlabel(var2)
    plt.plot(y[baseRec1.best_index], z[baseRec1.best_index], marker='o', markersize=10,mfc='none', color='k', zorder=10)
    if n_iter_folders==1:
        plt.axvline(x=std_value2,label='standard_value', linestyle='--',alpha=0.5)
        plt.plot(y[baseRec1.best_index], z[baseRec1.best_index], marker='o', markersize=10,mfc='none', color='k', label='best values',zorder=10)
    plt.legend( bbox_to_anchor=(leg_x, leg_y ) )
   
    ax_all3=plt.subplot(223)
    ax_all3.set_title("mod. factor phi1 vs "+var1)
    plt.xlabel(var1)
    plt.errorbar(x,z1,yerr=z1_err, fmt='o',label=LW_folder )
    plt.plot(x[baseRec1.best_index], z1[baseRec1.best_index], marker='o',mfc='none', markersize=10, color='k',zorder=10)
    if n_iter_folders==1:
          plt.axvline(x=std_value1,label='standard_value', linestyle='--',alpha=0.5)
          plt.plot(x[baseRec1.best_index], z1[baseRec1.best_index], marker='o', markersize=10,mfc='none',   color='k', label='best values',zorder=10)
    
    plt.legend( bbox_to_anchor=(leg_x, leg_y ) )

    ax_all4=plt.subplot(224)
    ax_all4.set_title("mod. factor phi1 vs "+var2)
    plt.errorbar(y,z1,yerr=z1_err, fmt='o',label=LW_folder )
    plt.xlabel(var2)
    plt.plot(y[baseRec1.best_index], z1[baseRec1.best_index], marker='o',mfc='none',  markersize=10, color='k', zorder=10)
    if n_iter_folders==1:
        plt.axvline(x=std_value2,label='standard_value', linestyle='--',alpha=0.5)
        plt.plot(y[baseRec1.best_index], z1[baseRec1.best_index], marker='o',mfc='none', markersize=10, color='k', label='best values',zorder=10)
    
    plt.legend( bbox_to_anchor=(leg_x, leg_y ) )
   # plt.show()


   

    
    outfilePlot=base_dir+'summary3.png'
    print ("outFile png =",outfilePlot)
    figAll_phi2.savefig(outfilePlot)

  
    


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

std_val1=-100
std_val2=-100


n_raw_opt=[]
n_ecut_opt=[]

n_raw_std=[]
n_ecut_std=[]

resolution_opt=[]
resolution_std=[]


n_iter_folders=0

#inizio scan su zero_thr

for LW_index in range (0,LW_bins):
#for LW_index in range (0,1):
        

    LW_folder='LWbin_'+str(LW_index)

    n_iter_folders+=1
    out_dir=base_dir+'/scan_results/'+LW_folder+'/'
    run_command('mkdir -p '+out_dir)
    
    baseRec1=base_rec()
    
    for  i in range (1,n_iters+1):
        #work_dir=out_dir+str(i)
        work_dir=base_dir+'/rec_'+str(i)+'/'+LW_folder
        file_out=work_dir+'/prova_outLW.txt'
        file_cfg= base_dir+'/rec_'+str(i) +'/config_simo.txt'
        
        print('reading file:',file_out)  
        baseRec1.read_file_rec(file_out,file_cfg)

  
    










     
    ##############################



        
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
    
    
    
    
    #title='energy= '+str(dict_energy[folder])+' KeV  (folder '+folder+')' 
    title=LW_folder 
    plot_all(baseRec1,out_dir,LW_folder) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    #fill variables for final plots (vs energy)
    #energy.append(dict_energy[folder])
    energy.append(baseRec1.LWmean[0]   ) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    best_var1.append (baseRec1.dict_rec[var1][baseRec1.best_index] )
    best_var1_up.append (baseRec1.max_x )
    best_var1_low.append (baseRec1.min_x )

    best_var2.append (baseRec1.dict_rec[var2][baseRec1.best_index] )
    best_var2_up.append (baseRec1.max_y )
    best_var2_low.append (baseRec1.min_y )

    #std_var1.append (baseRec1.dict_rec[var1][baseRec1.best_index] )
    
    
    n_raw_opt.append(baseRec1.dict_rec['n_raw'][baseRec1.best_index])
    #n_ecut_opt.append(baseRec1.dict_rec['n_ecut'][baseRec1.best_index])
    mod1.append(baseRec1.dict_rec['modulation1'][baseRec1.best_index])
    mod1_err.append(baseRec1.dict_rec['modulation1_err'][baseRec1.best_index])
    mod2.append(baseRec1.dict_rec['modulation2'][baseRec1.best_index])
    mod2_err.append(baseRec1.dict_rec['modulation2_err'][baseRec1.best_index])
   # resolution_opt.append(baseRec1.dict_rec['resolution2'][baseRec1.best_index])

    

    mod1std.append(baseRec1.dict_rec['modulation1'][std_index])
    mod1std_err.append(baseRec1.dict_rec['modulation1_err'][std_index])
    mod2std.append(baseRec1.dict_rec['modulation2'][std_index])
    mod2std_err.append(baseRec1.dict_rec['modulation2_err'][std_index])
    #resolution_std.append(baseRec1.dict_rec['resolution2'][std_index])
    n_raw_std.append(baseRec1.dict_rec['n_raw'][std_index])
   # n_ecut_std.append(baseRec1.dict_rec['n_ecut'][std_index])

    std_val1=baseRec1.dict_rec[var1][std_index]
    std_val2=baseRec1.dict_rec[var2][std_index]
    


###############
# write outFile.txt


outFileName=base_dir+'outScan_'+var1+'-'+var2+'.txt'
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

outFile.write( ' '.join(map(str,best_var1))+'\n'  )
outFile.write( ' '.join(map(str,best_var1_up))+'\n'  )
outFile.write( ' '.join(map(str,best_var1_low))+'\n'  )

outFile.write( ' '.join(map(str,best_var2))+'\n'  )
outFile.write( ' '.join(map(str,best_var2_up))+'\n'  )
outFile.write( ' '.join(map(str,best_var2_low))+'\n'  )

outFile.close()






    


# final plots:

##########
#  FIG 1
##########
"""
fig1=plt.figure(figsize=(20,10))
fig1.suptitle("n. events", fontsize=16)
fig1.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

# PLOT  mod_2 vs zero threshold
ax01=plt.subplot(221)
ax01.set_title('Event raw number')
plt.errorbar(energy,  n_raw_opt ,  fmt='bo--',label='n_raw ')
plt.xlabel('energy [KeV]')
plt.ylabel('n. raw events')
plt.legend()

ax02=plt.subplot(222)
ax02.set_title('n. events Ecuts')
plt.errorbar(energy,  n_ecut_opt,  fmt='bo--',label='n_ecut opt ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio n. ecuts')
plt.legend()



ax03=plt.subplot(223)
ax03.set_title('n. Ecuts/n_raw')
plt.errorbar(energy,  np.array(n_ecut_opt)/np.array(n_raw_opt),  fmt='bo--',label='n_ecut/n_raw opt ')
plt.errorbar(energy, np.array(n_ecut_std)/np.array(n_raw_std),  fmt='ro--',label='n_ecut/n_raw  std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio n. events')
plt.legend()

ax04=plt.subplot(224)
ax04.set_title('ratio  n. Ecuts/n_raw   opt/std ')
plt.errorbar(energy,  ( np.array(n_ecut_opt)/np.array(n_raw_opt))/(np.array(n_ecut_std)/np.array(n_raw_std))    ,  fmt='bo--',label='event fraction  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()



outfilePlot=base_dir+'summary1.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

"""



#############
#  FIG 2
##############

fig2=plt.figure(figsize=(20,10))
fig2.suptitle(var1+" - "+var2+" scan", fontsize=16)
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


# rapporti modulazione
ax2=plt.subplot(234)

#creo numpy array dalle liste
mod1_np=np.array(mod1)
mod2_np=np.array(mod2)
mod1Err_np=np.array(mod1_err)
mod2Err_np=np.array(mod2_err)

mod1std_np=np.array(mod1std)
mod2std_np=np.array(mod2std)
mod1stdErr_np=np.array(mod1std_err)
mod2stdErr_np=np.array(mod2std_err)


ratio1=mod1_np/mod1std_np
ratio1Err=(  (1./mod1std_np**2)*(mod1Err_np**2)+((mod1_np/(mod1std_np**2))**2)*(mod1stdErr_np**2) )**0.5

ratio2=mod2_np/mod2std_np
ratio2Err=(  (1./mod2std_np**2)*(mod2Err_np**2)+((mod2_np/(mod2std_np**2))**2)*(mod2stdErr_np**2) )**0.5


plt.errorbar(energy, ratio1, yerr=ratio1Err,  fmt='bo--',label='phi1  ratio opt/std ')
plt.errorbar(energy, ratio2, yerr=ratio2Err,  fmt='ro--',label='phi2  ratio opt/std ')
plt.xlabel('energy [KeV]')
plt.ylabel('ratio')
plt.legend()

#plots resolution

#ax3=plt.subplot(232)
#ax3.set_title('Resolution')
#plt.errorbar(energy, resolution_opt,  fmt='bo--',label='resolution opt ')
#plt.errorbar(energy, resolution_std,  fmt='ro--',label='resolution std ')
#plt.xlabel('energy [KeV]')
#plt.ylabel('resolution')
#plt.legend()

#plt.subplot(235)
#plt.errorbar(energy, np.array(resolution_opt)/np.array(resolution_std),  fmt='bo--',label='resolution  ratio opt/std ')
#plt.xlabel('energy [KeV]')
#plt.ylabel('ratio')
#plt.legend()

ax4=plt.subplot(233)
#best var1
ax4.set_title('best '+var1)
plt.errorbar(energy,best_var1, fmt='bo',label='best value')
ax4.fill_between(energy,best_var1_low, best_var1_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax4.axhline(y=std_val1,label='standard_value', linestyle='--',alpha=0.5)
plt.xlabel('energy [KeV]')
plt.ylabel('best '+var1)
plt.legend()
#plt.rc('grid',axes=True,  linestyle=":", color='grey')
#plt.grid(True,linestyle=':', color='grey')

# best var2
ax5=plt.subplot(236)
ax5.set_title('best '+var2)
#plt.errorbar(energy,n_raw_opt, fmt='bo--')
plt.errorbar(energy,best_var2, fmt='bo',label='best value')
ax5.fill_between(energy,best_var2_low, best_var2_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax5.axhline(y=std_val2,label='standard_value', linestyle='--',alpha=0.5)
plt.xlabel('energy [KeV]')
plt.ylabel('best '+var2)
plt.legend()



outfilePlot=base_dir+'summary2.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

                 
plt.show()


#plt.show()
 

