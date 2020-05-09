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






######################

base_dirPolarized='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/'



def plotDeltas(fileSummaryPol, fileSummarySpurius, parameter, yLims,  outFile, x_scPol,x_scUnPol,deltaPol_sc,deltaSpuria_sc):

    usa_12=0  # flag to use phi1 at E<3  (=1 => usa phi1)

    data1Moma1_pol=readSummaryData_1d(fileSummaryPol)
    data1Moma1_unpolDff=readSummaryData_unpolDFF_1d(fileSummarySpurius)
    
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax=plt.subplot(111)
     
    plt.grid(True,linestyle=':', color='grey') 

    ax.set_title('mod. differences '+ parameter)




 
    polStd12=np.ndarray(len (data1Moma1_pol.energy))
    
    polStd12Err=np.ndarray(len (data1Moma1_pol.energy))

    spuriaStd12=np.ndarray(len ( data1Moma1_unpolDff.energy))
    spuriaStd12Err=np.ndarray(len ( data1Moma1_unpolDff.energy))
 
    for i in range(0, len (data1Moma1_pol.energy)):
                   if data1Moma1_pol.energy[i]>3 :  # se E>3 uso phi2      
                      polStd12[i]= data1Moma1_pol.mod2std[i]
                      polStd12Err[i]= data1Moma1_pol.mod2std_err[i]
                      print("i=",i," val =",  data1Moma1_pol.mod2std[i], "polstd12[i]=",  polStd12[i] )                   
                   else:
                       polStd12[i]= data1Moma1_pol.mod1std[i]
                       polStd12Err[i]= data1Moma1_pol.mod1std_err[i]
                      
    for i in range(0, len (data1Moma1_unpolDff.energy)):
                   if data1Moma1_unpolDff.energy[i]>3 :  # se E>3 uso phi2      
                      spuriaStd12[i]=data1Moma1_unpolDff.mod2std[i]
                      spuriaStd12Err[i]=data1Moma1_unpolDff.mod2std_err[i]
                                        
                   else:
                       spuriaStd12[i]=data1Moma1_unpolDff.mod1std[i]
                       spuriaStd12Err[i]=data1Moma1_unpolDff.mod1std_err[i]
    

    #debug:
    print ("data1Moma1_pol.energy=",data1Moma1_pol.energy)
    print ("data1Moma1_pol.mod2std=",data1Moma1_pol.mod2std)
    print ("data1Moma1_pol.mod1std=",data1Moma1_pol.mod1std)
    print ("polStd12=",polStd12)
     
                       
                   
    delta_pol=(data1Moma1_pol.mod2-data1Moma1_pol.mod2std)
    delta_polErr=(data1Moma1_pol.mod2_err**2+data1Moma1_pol.mod2std_err**2)**0.5

    delta_polSpuria=(data1Moma1_unpolDff.mod_phi2_bestPol - data1Moma1_unpolDff.mod2std)
    delta_polSpuriaErr=(data1Moma1_unpolDff.mod_phi2_bestPolErr**2 + data1Moma1_unpolDff.mod2std_err**2)**0.5



    if usa_12==1:   
        delta_pol=(data1Moma1_pol.mod2-polStd12)
        delta_polErr=(data1Moma1_pol.mod2_err**2+polStd12Err**2)**0.5

        delta_polSpuria=(data1Moma1_unpolDff.mod_phi2_bestPol - spuriaStd12)
        delta_polSpuriaErr=(data1Moma1_unpolDff.mod_phi2_bestPolErr**2 + spuriaStd12Err**2)**0.5

    

                   
    plt.errorbar(data1Moma1_pol.energy, delta_pol, yerr=delta_polErr, fmt='ro',label='polarized data')

   
    plt.errorbar(data1Moma1_unpolDff.energy, delta_polSpuria, yerr=delta_polSpuriaErr, fmt='bo',label='spurious mod.')

    if yLims[0]!=yLims[1]:   # set limits only of the two parameters are different
        ax.set_ylim(yLims[0],yLims[1])

    plt.xlabel('energy [KeV]')
    plt.ylabel('Mod_bestPar - Mod_stdPar')

    # punti da slieds Simone Castellano (ove disponibili)
    if len(x_scPol)>0:
         plt.plot(x_scPol, deltaPol_sc, 'ro',mfc='none',markersize=10, label='polarized  sc')
         plt.plot(x_scUnPol, deltaSpuria_sc, 'bo',mfc='none',markersize=10, label='spurious mod.  sc')
         


    
    plt.legend()
    reorderLegend(ax)
    outfilePlot=outFile
    print ("outFile png =",outfilePlot)
    fig.savefig(outfilePlot)








##########################

#fileSummaryPol='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/


fileSummarySpuriusBest='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/bestPara/outBestPara.txt'

fileSummarySpurius1d='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanDmin/outScan_dmin.txt'

#outPlotFile=base_dirPolarized+'/scanMoma1/diff_Moma1.png'




data_unpolDff1d=readSummaryData_unpolDFF_1d(fileSummarySpurius1d)

data_unpolDffBest=readSummaryData_1d(fileSummarySpuriusBest)



x=data_unpolDffBest.energy


spuriaBest=(data_unpolDffBest.mod2)
spuriaBest_err=(data_unpolDffBest.mod2_err)


spuriaStd1=(data_unpolDff1d.mod1std)
spuriaStd1_err=(data_unpolDff1d.mod1std_err)

spuriaStd2=(data_unpolDff1d.mod2std)
spuriaStd2_err=(data_unpolDff1d.mod2std_err)




delta1=spuriaBest-spuriaStd1
delta1_err=(spuriaBest_err**2+spuriaStd1_err**2)**0.5

delta2=spuriaBest-spuriaStd2
delta2_err=(spuriaBest_err**2+spuriaStd2_err**2)**0.5


 
fig=plt.figure(figsize=(10,7))
fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
ax=plt.subplot(111)
     
plt.grid(True,linestyle=':', color='grey') 

ax.set_title('spuruius modulation: best - standard')


plt.errorbar(x, delta1, yerr=delta1_err, fmt='ro',label='phi 1')
plt.errorbar(x, delta2, yerr=delta2_err, fmt='bo',label='phi 2')

plt.xlabel('energy [KeV]')
plt.ylabel('Mod_bestPar - Mod_stdPar')

plt.legend()
reorderLegend(ax)
outfilePlot='diff_spuriaBest-std.png'
print ("outFile png =",outfilePlot)
fig.savefig(outfilePlot)




fig2=plt.figure(figsize=(10,7))
fig2.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
ax2=plt.subplot(111)
     
plt.grid(True,linestyle=':', color='grey') 

ax2.set_title('spurious mod. - best/std')

ratio1=spuriaBest/spuriaStd1
ratio2=spuriaBest/spuriaStd2
plt.errorbar(x, ratio1,  fmt='ro',label='phi 1')
plt.errorbar(x, ratio2,  fmt='bo',label='phi 2')

plt.xlabel('energy [KeV]')
plt.ylabel('Mod_bestPar - Mod_stdPar')

plt.legend()
reorderLegend(ax)
outfilePlot='diff_spuriaBest-std.png'
print ("outFile png =",outfilePlot)
fig.savefig(outfilePlot)






plt.show()


