# loop su scan ixperecon 1dim (variando un solo parametro) - versione per MODULAZIONE SPURIA 
# produce i plot a energia fissa ( modulazione, risoluzione, chi2, etc al variare del parametro)
# produce i plot riassuntivi finali (parametro ottimale, modulazione best, etc) vs energia
# produce file riassuntivo di output



from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl



#import plotSummary.plotData1d



#import matplotlib.cm as cm
#from matplotlib.colors import Normalize



mpl.rcParams['legend.loc'] = 'upper right'   # default position
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True
mpl.rcParams['font.size']=15  #!!!!!!!!!!!!!!!!!!!!!!!!!!



# PARAMETERS:

first_iter=1

""""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanZeroThr/'
x_var='zero_thr'
std_index=9-1
n_iters=19
maximize='both'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/outScan_zero_thr.txt'
loc_ratios='upper right'  # posizione legenda plot rapporti
loc_deltas='upper right'  #  "          "      "   differenze
"""



"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanMoma1/'
x_var='moma1_thr'
std_index=9-1
n_iters=20
maximize='both'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/outScan_moma1_thr.txt'
loc_ratios='upper right'
loc_deltas='upper right'
"""


"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanMoma2/'
x_var='moma2_thr'
std_index=9-1
n_iters=20
maximize='phi2'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/outScan_moma2_thr.txt'
loc_ratios='upper right'
loc_deltas='upper right'
"""


"""
base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanDmin/'
x_var='dmin'
std_index=8-1
n_iters=17
maximize='phi2'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin/outScan_dmin.txt'
loc_ratios='upper left'
loc_deltas='lower left'
"""



base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scan_ws/'
x_var='weight_scale'
std_index=5-1-1
n_iters=20
maximize='phi2'  # phi1, ph2, both
summaryFile='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/outScan_weight_scale.txt'
loc_ratios='upper right'
loc_deltas='lower right'
first_iter=2




energy_dirs=['2.04','2.29', '2.70', '2.98', '3.69', '5.89']
#energy_dirs=['2.04','2.29', '2.70', '2.98']

#energy_dirs=['5.89']
eps_dir={'2.04':['000647','000658']  ,'2.29':['000669','000677'], '2.70':['000686','000691'], '2.98':['000704', '000726'], '3.69':['000733','000740'], '5.89':['000744','000752']}   # dir angolo1 e angolo2





class plotData1d:


   def __init__(self,fileName): 
       self.fileName=fileName

       self.readFile()


   def readFile(self):   
       infile=open(self.fileName,'r')
       npArrayList=[]

       for line in infile:
           #print(line[:-1])
           splittedLine=line[:-1].split(" ")
           #floatList=list(map(float, splittedLine))
           #print ("flaotList = " ,floatList )
           #floatList2=[float(i) for i in splittedLine ]
           #print ("float2=",floatList2 )

           npArrayList.append(np.array(splittedLine,dtype=np.float32))
           #print ("npArray=",npArray )


       self.energy=npArrayList[0]
       self.mod1=npArrayList[1]
       self.mod1_err=npArrayList[2]
       self.mod2=npArrayList[3]
       self.mod2_err=npArrayList[4]
       self.mod1std=npArrayList[5]
       self.mod1std_err=npArrayList[6]
       self.mod2std=npArrayList[7]
       self.mod2std_err=npArrayList[8]
       self.best_val=npArrayList[9]
       self.best_val_up=npArrayList[10]
       self.best_val_low=npArrayList[11]


       print ("energy =",self.energy)
       print ("best_val_low=",self.best_val_low)
       print ("best_val=",self.best_val)
















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




def plot_all(baseRec1,outdir,folder):

    title='energy= '+folder+' KeV' 
    
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
    figAll_phi2=plt.figure(1,figsize=(20,10) )
    figAll_phi2.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all=plt.subplot(121)
    ax_all.set_title("mod. factor phi2")
    plt.xlabel(x_var)
    plt.errorbar(x,y,yerr=y_err, fmt='o--',label='E='+folder+' KeV' )
    plt.legend( bbox_to_anchor=(1.1, 1.11) )

    ax_all2=plt.subplot(122)
    ax_all2.set_title("mod. factor phi1")
    plt.errorbar(x,y2,yerr=y2_err, fmt='o--',label='E='+folder+' KeV' )
    plt.xlabel(x_var)
    plt.legend(bbox_to_anchor=(1.1, 1.11))


    outfilePlot=base_dir+'summary3.png'
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
    
    
    # plot cumulativo (all E), con mod relativa ai valori nominali  PHI2

    figAll_rel=plt.figure(2,figsize=(10,7) )
    figAll_rel.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    ax_all=plt.subplot(111)
    ax_all.set_title("mod. factor phi2 / mod factor phi2 [standard_parameters] ")
    plt.xlabel(x_var)

    ratio=mod1_np/nom_mod
    #ratioErr=(  (1./nom_mod**2)*(mod1Err_np**2)+((mod1_np/(nom_mod**2))**2)*(nom_mod_err**2) )**0.5
    #plt.errorbar(x,ratio,yerr=ratioErr, fmt='o--',label='E='+folder+' KeV' )
    plt.errorbar(x,ratio, fmt='o--',label='E='+folder+' KeV' )
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
    plt.errorbar(x,ratio2, fmt='o--',label='E='+folder+' KeV' )
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
    #ratioErr=(  (1./nom_mod**2)*(mod1Err_np**2)+((mod1_np/(nom_mod**2))**2)*(nom_mod_err**2) )**0.5
    #plt.errorbar(x,ratio,yerr=ratioErr, fmt='o--',label='E='+folder+' KeV' )
    plt.errorbar(x,delta, fmt='o--',label='E='+folder+' KeV' )
    #plt.legend( bbox_to_anchor=(1.1, 1.11) )
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

    #mod2_np=np.array(y2)
    #mod2Err_np=np.array(y2_err)
    #nom_mod2=mod2_np[std_index]
    #nom_mod2_err= mod2Err_np[std_index]

    delta2=mod2_np-nom_mod2
    #ratioErr=(  (1./nom_mod**2)*(mod1Err_np**2)+((mod1_np/(nom_mod**2))**2)*(nom_mod_err**2) )**0.5
    #plt.errorbar(x,ratio,yerr=ratioErr, fmt='o--',label='E='+folder+' KeV' )
    plt.errorbar(x,delta2, fmt='o--',label='E='+folder+' KeV' )
    #plt.legend( bbox_to_anchor=(1.1, 1.11) )
    plt.legend(loc=loc_deltas)
   
    outfilePlot=base_dir+'summary_deltas_phi1.png'
    print ("outFile png =",outfilePlot)
    figAll_delta2.savefig(outfilePlot)

    
    

    ###########################
    # plot singole energie

    ### prendo punto da summary
    summaryData1d=plotData1d(summaryFile)
    index_summary=1e6
    try:
       index_summary=np.where(  np.logical_and( summaryData1d.energy<(float(folder)+0.05),  summaryData1d.energy >(float(folder)-0.05) ) )[0][0]
    except:
       print ("print energia non trovata nello scan ploarizzato")
    
    fig01=plt.figure(figsize=(10,7))
    fig01.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
    fig01.suptitle('scan '+x_var, fontsize=16)
    # PLOT  mod_2 vs zero threshold

    ax01=plt.subplot(111)
    
    #ax01.set_title('scan '+x_var)
    plt.errorbar(x,y,yerr=y_err, fmt='bo-',label="modulation phi2")
    plt.errorbar(x,y2,yerr=y2_err, fmt='ro-',label="modulation phi1", alpha=0.3) # modulazione phi1
    
    plt.xlabel(x_var)
    plt.ylabel('modulation')

    plt.axvline(x=std_value,label='standard_value', linestyle='--')
    if (index_summary<1000):
       plt.axvline(x=summaryData1d.best_val[index_summary],label='best_value (data) ', linestyle='--', color='red')

    plt.legend(loc='upper right')
    outfile=out_dir+'mod_both.png'
    print ("outFile png =",outfile)
    plt.savefig(outfile)



    
    # PLOT  mod_1 vs zero threshold
    #fig02=plt.figure(figsize=(10,7))
    #fig01.suptitle(title, fontsize=16)
    #fig02.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)

    
    #ax02=plt.subplot(111)
    #ax02.set_title('modulation phi_1')
    #plt.errorbar(x,y2,yerr=y2_err, fmt='ro--',label="modulation phi1") # modulazione phi1
    #plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.9)
    #plt.legend(loc='lower right')
    #plt.xlabel(x_var)
    #plt.ylabel('modulation_phi1')
    #if (index_summary<1000):
    #   plt.axvline(x=summaryData1d.best_val[index_summary],label='best_value (data) ', linestyle='--',alpha=0.9)
    #outfile=out_dir+'mod_phi1.png'
    #print ("outFile png =",outfile)
    #plt.savefig(outfile)




#   summary totale:
    fig03=plt.figure(figsize=(20,10))
    #fig01.suptitle(title, fontsize=16)
    fig03.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.09,hspace=0.250)
   
    
    # PLOT  mod_2 vs zero threshold
    ax07=plt.subplot(231)
    ax07.set_title('modulation phi_2')
    plt.errorbar(x,y,yerr=y_err, fmt='bo--',label="modulation phi2")
    plt.xlabel(x_var)
    plt.ylabel('modulation_2')

    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.9)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )

    if (index_summary<1000):
       plt.axvline(x=summaryData1d.best_val[index_summary],label='best_value (data) ', linestyle='--',alpha=0.9, color='red')
    if baseRec1.best_phi==2:
        ax07.add_patch(rect)
    plt.legend(loc='lower right')

    # mod phi1:
    ax08=plt.subplot(232)
    ax08.set_title('modulation phi_1')
    plt.errorbar(x,y2,yerr=y2_err, fmt='ro--',label="modulation phi1") # modulazione phi1
    plt.axvline(x=std_value,label='standard_value', linestyle='--',alpha=0.5)
    plt.axvline(x=best_thr,label='best value', linestyle='--',color='red',alpha=0.5 )
    plt.axvline(x=min_thr, linestyle=':',color='grey',alpha=0.5,label='band' )
    plt.axvline(x=max_thr, linestyle=':',color='grey',alpha=0.5 )
    if baseRec1.best_phi==1:
        ax08.add_patch(rect)
    plt.legend(loc='lower right')
    plt.xlabel(x_var)
    plt.ylabel('modulation_phi1')
    if (index_summary<1000):
       plt.axvline(x=summaryData1d.best_val[index_summary],label='best_value (data) ', linestyle='--',alpha=0.5)
  
    # plot chi2  
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

    # resolution
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

    # rapporto 'ratio n_physical/n_raw'
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

    # rapproto n_ecut/n_physical
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
    
    for  i in range (first_iter,n_iters+1): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    #for  i in range (2,n_iters+1):
       
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


    
    title='energy= '+folder+' KeV' 
    plot_all(baseRec1,out_dir,folder) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
print('test string=',' '.join(map(str,energy)) )

###############
# write outFile.txt

outFileName=base_dir+'outScan_'+x_var+'.txt'
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


    
outfilePlot=base_dir+'summary1.png'
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
plt.legend()


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




outfilePlot=base_dir+'summary2.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

                 
plt.show()
