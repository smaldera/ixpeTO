from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

from readSummaryData import *


# spostare da qualche parte....


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

base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/'



data1Dmin=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin/outScan_dmin_v2.txt')
data1ws=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/outScan_weight_scale_v2.txt')

data2wsDmin=plotData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/outScan_dmin-weight_scale.txt')


#tutti i plot per rapporto 
data1ZeroThr=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/outScan_zero_thr_v2.txt')
data1Moma1=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/outScan_moma1_thr_v2.txt')
data1Moma2=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/outScan_moma2_thr_v2.txt')
data1Dmax=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmax/outScan_dmax_v2.txt')

data2zeroMoma=plotData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroMoma12/outScan_zero_thr-moma1_thr.txt')


######################
# print obtimal values:

print("\n \n BEST PARAMETERS")
print("bestZeroThr=",data2zeroMoma.best_var1)
print("bestMoma12=",data2zeroMoma.best_var2)
print("bestDmin=",data2wsDmin.best_var1)
print("bestws=",data2wsDmin.best_var2)


print("E=",data2wsDmin.energy)





############################3
# PLOT: BEST Dmin vs energy

# dati s. castellano (from scan on MC sims):
eMC = np.array( [2.01, 2.29, 2.70, 2.98, 3.69, 4.00, 4.50, 5.00, 6.00, 6.40, 7.00, 8.00])## energy keV
dminMC = np.array([0.3, 0.3, 0.7, 0.7, 1.7, 2.0, 2.1, 2.1, 1.9, 1.9, 1.8, 1.7])## dmin best




plt.figure(1)                     
ax1=plt.subplot(111)
ax1.set_title('best dmin')
plt.errorbar(data1Dmin.energy,data1Dmin.best_val, fmt='bo',label='best dmin')
#ax1.fill_between(data1Dmin.energy,data1Dmin.best_val_low, data1Dmin.best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax1.fill_between(data1Dmin.energy,data1Dmin.best_val_low, data1Dmin.best_val_up,color='gray',alpha=0.1, interpolate=True)

ax1.axhline(y=1.5,label='standard_value', linestyle='--',alpha=0.5)
plt.errorbar(data2wsDmin.energy,data2wsDmin.best_var1, fmt='ro',label='best dmin 2d scan')
#ax1.fill_between(data2wsDmin.energy,data2wsDmin.best_var1_low, data2wsDmin.best_var1_up,color='red',alpha=0.1, interpolate=True,label=r'1$\sigma$ band 2d scan')
ax1.fill_between(data2wsDmin.energy,data2wsDmin.best_var1_low, data2wsDmin.best_var1_up,color='red',alpha=0.1, interpolate=True)


plt.plot(eMC,dminMC,'ko',mfc='none' , label='best dmin MC')

plt.xlabel('energy [KeV]')
plt.ylabel('best dmin')
plt.legend()
#reorderLegend(ax1,['best dmin',r'1$\sigma$ band' ,'best dmin 2d scan', r'1$\sigma$ band 2d scan',  'best dmin MC', 'standard_value' ])
reorderLegend(ax1,['best dmin','best dmin 2d scan',  'best dmin MC', 'standard_value' ])



outfilePlot=base_dir+'Dmin_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

###################
# fig2# PLOT: BEST weight scale vs energy

# dati MC (s.castellano)
eMC_ws = np.array([2.01, 2.29, 2.70, 2.98, 3.00, 3.69, 4.50, 6.40])## energy keV
wsMC = np.array([0.16, 0.20, 0.14, 0.19, 0.14, 0.04, 0.03, 0.03])#ws best -> weight scale


plt.figure(2)                     
ax2=plt.subplot(111)
ax2.set_title('best weight scale')
plt.errorbar(data1ws.energy,data1ws.best_val, fmt='bo',label='best ws')
#ax2.fill_between(data1ws.energy,data1ws.best_val_low, data1ws.best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax2.fill_between(data1ws.energy,data1ws.best_val_low, data1ws.best_val_up,color='gray',alpha=0.1, interpolate=True)

ax2.axhline(y=0.05,label='standard_value', linestyle='--',alpha=0.5)


plt.errorbar(data2wsDmin.energy,data2wsDmin.best_var2, fmt='ro',label='best ws 2d scan')
#ax2.fill_between(data2wsDmin.energy,data2wsDmin.best_var2_low, data2wsDmin.best_var2_up,color='red',alpha=0.1, interpolate=True,label=r'1$\sigma$ band 2d scan'   )
ax2.fill_between(data2wsDmin.energy,data2wsDmin.best_var2_low, data2wsDmin.best_var2_up,color='red',alpha=0.1, interpolate=True   )


plt.plot(eMC_ws,wsMC,'ko',mfc='none',markersize=10, label='best ws MC')
plt.xlabel('energy [KeV]')
plt.ylabel('best ws')
plt.legend()
reorderLegend(ax2,['best ws','best ws 2d scan',  'best ws MC', 'standard_value' ])

outfilePlot=base_dir+'ws_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

##################################################33
# fig 4
# PLOT Dmax

# dati scan MC (s. castellano)
eMC_dmax = np.array([2.01, 2.29, 2.70, 2.98, 3.69, 4.50, 6.40])## energy keV
dmaxMC =np.array( [2.3, 2.3, 2.3, 2.3, 3.9, 4.7, 4.3])#dmax best

plt.figure(4)                     
ax4=plt.subplot(111)
ax4.set_title('best dmax')
plt.errorbar(data1Dmax.energy,data1Dmax.best_val, fmt='bo',label='best Dmax')
#ax4.fill_between(data1Dmax.energy,data1Dmax.best_val_low, data1Dmax.best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax4.fill_between(data1Dmax.energy,data1Dmax.best_val_low, data1Dmax.best_val_up,color='gray',alpha=0.1, interpolate=True)

ax4.axhline(y=3.5,label='standard_value', linestyle='--',alpha=0.5)

plt.plot(eMC_dmax,dmaxMC,'ko',mfc='none',markersize=10, label='best Dmax MC')
plt.xlabel('energy [KeV]')
plt.ylabel('best Dmax')
plt.legend()
reorderLegend(ax4,['best Dmax', 'best Dmax MC', 'standard_value' ])



outfilePlot=base_dir+'dmax_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)


##############################################
# fig 5
# PLOT moma1

# dati scan MC(s.castellano)
moma1MC =np.array( [22,28,34,32,22,20,20])#moma1 best (energie tue senza aggiunte)




plt.figure(5)                     
ax5=plt.subplot(111)
ax5.set_title('best moma1')
plt.errorbar(data1Moma1.energy,data1Moma1.best_val, fmt='bo',label='best Moma1')
ax5.fill_between(data1Moma1.energy,data1Moma1.best_val_low, data1Moma1.best_val_up,color='gray',alpha=0.1, interpolate=True)
ax5.axhline(y=36,label='standard_value', linestyle='--',alpha=0.5)


plt.errorbar(data2zeroMoma.energy,data2zeroMoma.best_var2, fmt='ro',label='best moma1(2) 2D scan')
ax5.fill_between(data2zeroMoma.energy,data2zeroMoma.best_var2_low, data2zeroMoma.best_var2_up,color='red',alpha=0.1, interpolate=True )


plt.plot(data1Moma1.energy,moma1MC,'ko',mfc='none',markersize=10, label='best moma1 MC')

plt.xlabel('energy [KeV]')
plt.ylabel('best moma1')
plt.legend()
reorderLegend(ax5,['best Moma1', 'best moma1(2) 2D scan','best moma1 MC', 'standard_value' ])


outfilePlot=base_dir+'moma1_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

#########################################
# fig 6
# PLOT moma2

# dati scan MC(s.castellano)
moma2MC =np.array([30,20,24,20,22,20,20] )#moma1 best 

plt.figure(6)                     
ax6=plt.subplot(111)
ax6.set_title('best moma2')
plt.errorbar(data1Moma2.energy,data1Moma2.best_val, fmt='bo',label='best Moma2')
ax6.fill_between(data1Moma2.energy,data1Moma2.best_val_low, data1Moma2.best_val_up,color='gray',alpha=0.1, interpolate=True)
ax6.axhline(y=36,label='standard_value', linestyle='--',alpha=0.5)


plt.errorbar(data2zeroMoma.energy,data2zeroMoma.best_var2, fmt='ro',label='best moma1(2) 2D scan')
ax6.fill_between(data2zeroMoma.energy,data2zeroMoma.best_var2_low, data2zeroMoma.best_var2_up,color='red',alpha=0.1, interpolate=True  )
plt.plot(data1Moma2.energy,moma2MC,'ko',mfc='none',markersize=10, label='best moma2 MC')

plt.xlabel('energy [KeV]')
plt.ylabel('best moma2')
plt.legend()
reorderLegend(ax6,['best Moma2', 'best moma1(2) 2D scan','best moma2 MC', 'standard_value' ])




outfilePlot=base_dir+'moma2_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)


#########################################
# fig 7
# PLOT clustering threshold


plt.figure(7)                     
ax7=plt.subplot(111)
ax7.set_title('best clustering threshold')
plt.errorbar(data1ZeroThr.energy,data1ZeroThr.best_val, fmt='bo',label='best clustering thr.')
ax7.fill_between(data1ZeroThr.energy,data1Moma2.best_val_low, data1Moma2.best_val_up,color='gray',alpha=0.1, interpolate=True)
ax7.axhline(y=20,label='standard_value', linestyle='--',alpha=0.5)


plt.errorbar(data2zeroMoma.energy,data2zeroMoma.best_var1, fmt='ro',label='best clustering thr 2d scan')
ax7.fill_between(data2zeroMoma.energy,data2zeroMoma.best_var1_low, data2zeroMoma.best_var1_up,color='red',alpha=0.1, interpolate=True )


plt.xlabel('energy [KeV]')
plt.ylabel('best clustering threshold')
plt.legend()
reorderLegend(ax7,['best clustering thr.', 'best clustering thr 2d scan', 'standard_value' ])


outfilePlot=base_dir+'clusterThr_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)


##########################
# fig 3
plt.figure(3)
ax3=plt.subplot(111)
plt.grid(True,linestyle=':', color='grey') 

ax3.set_title('modPhi2best/ max(modPhi1std,modPhi2std)')
plt.errorbar(data1ws.energy,data1ws.mod2/(np.maximum(data1ws.mod2std,data1ws.mod1std)), fmt='bo',label='scan ws')
plt.errorbar(data1Dmin.energy,data1Dmin.mod2/(np.maximum(data1Dmin.mod2std,data1Dmin.mod1std)), fmt='ro',label='scan dmin')
plt.errorbar(data2wsDmin.energy,data2wsDmin.mod2/(np.maximum(data2wsDmin.mod2std,data2wsDmin.mod1std)), fmt='ko',label='scan dmin-ws')



plt.errorbar(data1ZeroThr.energy, (np.maximum(data1ZeroThr.mod2, data1ZeroThr.mod1 ))/(np.maximum(data1ZeroThr.mod2std,data1ZeroThr.mod1std)), fmt='go',label='scan clustering thr (max mo1d mod2)')
plt.errorbar(data1Moma1.energy,(np.maximum(data1Moma1.mod2,data1Moma1.mod1  ))/(np.maximum(data1Moma1.mod2std,data1Moma1.mod1std)), fmt='bs',label='scan moma1 thr  (max mod1 mod2)')
plt.errorbar(data1Moma2.energy,data1Moma2.mod2/(np.maximum(data1Moma2.mod2std,data1Moma2.mod1std)), fmt='rs',label='scan moma2 thr')
plt.errorbar(data1Dmax.energy,data1Dmax.mod2/(np.maximum(data1Dmax.mod2std,data1Dmax.mod1std)), fmt='r*',label='scan dmax')

plt.errorbar(data2zeroMoma.energy,(np.maximum(data2zeroMoma.mod2,data2zeroMoma.mod1 ))/(np.maximum(data2zeroMoma.mod2std,data2zeroMoma.mod1std)), fmt='kP',label='scan zero-moma12  (max mod1 mod2)  ')


plt.xlabel('energy [KeV]')
plt.ylabel('modPhi2best/ max(modPhi1std,modPhi2std)')
plt.legend()
outfilePlot=base_dir+'modComparison.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

plt.show()
