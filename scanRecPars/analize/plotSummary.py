from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt




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


       
class plotData2d:


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

       self.best_var1=npArrayList[9]
       self.best_var1_up=npArrayList[10]
       self.best_var1_low=npArrayList[11]

       self.best_var2=npArrayList[12]
       self.best_var2_up=npArrayList[13]
       self.best_var2_low=npArrayList[14]



#############################################################


base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/'


data1Dmin=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin/outScan_dmin.txt')
data1ws=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/outScan_weight_scale.txt')

data2wsDmin=plotData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/outScan_dmin-weight_scale.txt')


#tutti i plot per rapporto 
data1ZeroThr=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/outScan_zero_thr.txt')
data1Moma1=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/outScan_moma1_thr.txt')
data1Moma2=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/outScan_moma2_thr.txt')
data1Dmax=plotData1d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmax/outScan_dmax.txt')

data2zeroMoma=plotData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroMoma12/outScan_zero_thr-moma1_thr.txt')






plt.figure(1)                     
ax1=plt.subplot(111)
ax1.set_title('best dmin')
plt.errorbar(data1Dmin.energy,data1Dmin.best_val, fmt='bo',label='best dmin')
ax1.fill_between(data1Dmin.energy,data1Dmin.best_val_low, data1Dmin.best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax1.axhline(y=1.5,label='standard_value', linestyle='--',alpha=0.5)

plt.errorbar(data2wsDmin.energy,data2wsDmin.best_var1, fmt='ro',label='best dmin 2d scan  ')
ax1.fill_between(data2wsDmin.energy,data2wsDmin.best_var1_low, data2wsDmin.best_var1_up,color='red',alpha=0.1, interpolate=True,label=r'1$\sigma$ band 2d scan')



plt.xlabel('energy [KeV]')
plt.ylabel('best dmin')
plt.legend()
outfilePlot=base_dir+'Dmin_scanSummary.png'
print ("outFile png =",outfilePlot)
plt.savefig(outfilePlot)

###################
# fig2 
plt.figure(2)                     
ax2=plt.subplot(111)
ax2.set_title('best ws')
plt.errorbar(data1ws.energy,data1ws.best_val, fmt='bo',label='best ws')
ax2.fill_between(data1ws.energy,data1ws.best_val_low, data1ws.best_val_up,color='gray',alpha=0.1, interpolate=True,label=r'1$\sigma$ band')
ax2.axhline(y=0.05,label='standard_value', linestyle='--',alpha=0.5)


plt.errorbar(data2wsDmin.energy,data2wsDmin.best_var2, fmt='ro',label='best ws 2d scan  ')
ax2.fill_between(data2wsDmin.energy,data2wsDmin.best_var2_low, data2wsDmin.best_var2_up,color='red',alpha=0.1, interpolate=True,label=r'1$\sigma$ band 2d scan'   )


plt.xlabel('energy [KeV]')
plt.ylabel('best ws')
plt.legend()

outfilePlot=base_dir+'ws_scanSummary.png'
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
