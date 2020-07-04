from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from readSummaryData import *


# spostare da qualche parte....

mpl.rcParams['font.size']=15  #!!!!!!!!!!!!!!!!!!!!!!!!!!

#  Returns tuple of handles, labels for axis ax, after reordering them to conform to the label order `order`, and if unique is True, after removing entries with duplicate labels.
def reorderLegend(ax=None,order=None,unique=False):
    if ax is None: ax=plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0])) # sort both labels and handles by labels
    if order is not None: # Sort according to a given list (not necessarily complete)
        keys=dict(zip(order,range(len(order))))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
    if unique:  labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels)) # Keep only the first of each handle
    ax.legend(handles, labels)
    return(handles, labels)


def unique_everseen(seq, key=None):
    seen = set()
    seen_add = seen.add
    return [x for x,k in zip(seq,key) if not (k in seen or seen_add(k))]





def unisci_phi(energy, phi1,phi2):   # tutti np.array

    final=np.zeros(len(energy))
    
    for i in range (0,len(energy)):
        if energy[i]<3:
            final[i]=phi1[i]
        else:    
            final[i]=phi2[i]
            
    return final        

            
def plot_std(data_unpolDff1d,base_dir):


    x=data_unpolDff1d.energy
    spuriaStd1=(data_unpolDff1d.mod1std)
    spuriaStd1_err=(data_unpolDff1d.mod1std_err)

    spuriaStd2=(data_unpolDff1d.mod2std)
    spuriaStd2_err=(data_unpolDff1d.mod2std_err)

    spuriaStdFinal=unisci_phi(x,spuriaStd1,spuriaStd2 )
    spuriaStdFinal_err=unisci_phi(x,spuriaStd1_err,spuriaStd2_err )

    print("spuriaStdFinal=",spuriaStdFinal)
    

    # plot
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax=plt.subplot(111)
    
    plt.grid(True,linestyle=':', color='grey') 

    ax.set_title('spurious modulation: standard rec')

    plt.errorbar(x, spuriaStd1,yerr= spuriaStd1_err  , fmt='ro',label='phi 1')
    plt.errorbar(x,spuriaStd2, yerr=spuriaStd2_err, fmt='bo',label='phi 2')

    #plt.errorbar(x, spuriaStdFinal, 'ko', label='phi 1 E<=3KeV, phi2 E>3KeV')
    plt.errorbar(x, spuriaStdFinal, fmt='ko--',mfc='none',markersize=10 ,label='phi 1 E<=3KeV, phi2 E>3KeV')
   
  

    plt.xlabel('energy [KeV]')
    plt.ylabel('spurious modulation')

    plt.legend()
    #reorderLegend(ax)
    outfilePlot=base_dir+'/spuria_std.png'
    print ("outFile png =",outfilePlot)
    fig.savefig(outfilePlot)


    
def plot_best(data_unpolDffBest,base_dir):


    x=data_unpolDff1d.energy
    spuria1=(data_unpolDffBest.mod1std)
    spuria1_err=(data_unpolDffBest.mod1std_err)

    spuria2=(data_unpolDffBest.mod2std)
    spuria2_err=(data_unpolDffBest.mod2std_err)

    #spuriaStdFinal=unisci_phi(x,spuriaStd1,spuriaStd2 )
    #spuriaStdFinal_err=unisci_phi(x,spuriaStd1_err,spuriaStd2_err )
    
    
    # plot
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax=plt.subplot(111)
    
    plt.grid(True,linestyle=':', color='grey') 

    ax.set_title('spurious modulation: best rec (pha parametrization)  ')

    plt.errorbar(x, spuria1,yerr= spuria1_err  , fmt='ro',label='phi 1')
    plt.errorbar(x,spuria2, yerr=spuria2_err, fmt='bo--',label='phi 2')

    #plt.errorbar(x, spuriaStdFinal, 'ko', label='phi 1 E<=3KeV, phi2 E>3KeV')
    #plt.errorbar(x, spuriaStdFinal, fmt='ko',mfc='none',markersize=10 ,label='phi 1 E<=3KeV, phi2 E>3KeV')
   
  

    plt.xlabel('energy [KeV]')
    plt.ylabel('spurious modulation')

    plt.legend()
    #reorderLegend(ax)
    outfilePlot=base_dir+'/spuria_best.png'
    print ("outFile png =",outfilePlot)
    fig.savefig(outfilePlot)
    
       


    
def plot_ratios(data_unpolDff1d, data_unpolDffBest,base_dir):


    x=data_unpolDff1d.energy
    spuria1=(data_unpolDffBest.mod1std)
    spuria1_err=(data_unpolDffBest.mod1std_err)

    spuria2=(data_unpolDffBest.mod2std)
    spuria2_err=(data_unpolDffBest.mod2std_err)

  
    spuriaStd1=(data_unpolDff1d.mod1std)
    spuriaStd1_err=(data_unpolDff1d.mod1std_err)

    spuriaStd2=(data_unpolDff1d.mod2std)
    spuriaStd2_err=(data_unpolDff1d.mod2std_err)

    spuriaStdFinal=unisci_phi(x,spuriaStd1,spuriaStd2 )
    spuriaStdFinal_err=unisci_phi(x,spuriaStd1_err,spuriaStd2_err )



    
    
    # plot
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax=plt.subplot(111)
    
    plt.grid(True,linestyle=':', color='grey') 

    ax.set_title('Ratio best/std ')

    plt.errorbar(x, spuria2/spuriaStd1,  fmt='ro',label='phi2_best/ phi1 std')
    plt.errorbar(x,spuria2/spuriaStd2,  fmt='bo',label='phi2_best/ phi2 std')
    plt.errorbar(x, spuria2/spuriaStdFinal, fmt='ko--',mfc='none',markersize=10    ,label='phi2_best/(phi1_std @ E<3KeV,  ph2_std @ E>3KeV)')
    

    plt.xlabel('energy [KeV]')
    plt.ylabel('mod._best/mod_std')

    plt.legend()
    #reorderLegend(ax)
    outfilePlot=base_dir+'/Ratios_spuria.png'
    print ("outFile png =",outfilePlot)
    fig.savefig(outfilePlot)
    


def plot_diff(data_unpolDff1d, data_unpolDffBest,base_dir):


    x=data_unpolDff1d.energy
    spuria1=(data_unpolDffBest.mod1std)
    spuria1_err=(data_unpolDffBest.mod1std_err)

    spuria2=(data_unpolDffBest.mod2std)
    spuria2_err=(data_unpolDffBest.mod2std_err)

  
    spuriaStd1=(data_unpolDff1d.mod1std)
    spuriaStd1_err=(data_unpolDff1d.mod1std_err)

    spuriaStd2=(data_unpolDff1d.mod2std)
    spuriaStd2_err=(data_unpolDff1d.mod2std_err)

    spuriaStdFinal=unisci_phi(x,spuriaStd1,spuriaStd2 )
    spuriaStdFinal_err=unisci_phi(x,spuriaStd1_err,spuriaStd2_err )



    
    
    # plot
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax=plt.subplot(111)
    
    plt.grid(True,linestyle=':', color='grey') 

    ax.set_title('Difference best-std ')

    plt.errorbar(x, spuria2-spuriaStd1,  fmt='ro',label='phi2_best - phi1 std')
    plt.errorbar(x,spuria2-spuriaStd2,  fmt='bo',label='phi2_best - phi2 std')
    plt.errorbar(x, spuria2-spuriaStdFinal, fmt='ko--',mfc='none',markersize=10    ,label='phi2_best-(phi1_std@E<3KeV,  phi2_std@E>3KeV)')
    

    plt.xlabel('energy [KeV]')
    plt.ylabel('mod._best-mod_std')

    plt.legend()
    #reorderLegend(ax)
    outfilePlot=base_dir+'/Differences_spuria.png'
    print ("outFile png =",outfilePlot)
    fig.savefig(outfilePlot)
    




######################

base_dirPolarized='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/'




#fileSummaryPol='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/


fileSummarySpuriusBest='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/bestPara_phaParametrization/outBestPara_ixpereconParametrization.txt'
fileSummarySpurius1d='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanDmin/outScan_dmin.txt'



data_unpolDff1d=readSummaryData_unpolDFF_1d(fileSummarySpurius1d)
data_unpolDffBest=readSummaryData_1d(fileSummarySpuriusBest)



plot_std(data_unpolDff1d, base_dirPolarized) 
plot_best(data_unpolDffBest, base_dirPolarized) 

plot_ratios(data_unpolDff1d,  data_unpolDffBest, base_dirPolarized) 
plot_diff(data_unpolDff1d,  data_unpolDffBest, base_dirPolarized) 






plt.show()


