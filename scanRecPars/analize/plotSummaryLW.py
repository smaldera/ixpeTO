from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


from readSummaryData import *

mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['axes.grid'] = True
mpl.rcParams['font.size']=15  #!!!!!!!!!!!!!!!!!!!!!!!!!!


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



def unisci_phi(energy, phi1,phi2):   # tutti np.array

    final=np.zeros(len(energy))
    
    for i in range (0,len(energy)):
        if energy[i]<3:
            final[i]=phi1[i]
        else:    
            final[i]=phi2[i]
            
    return final        



def dmin_newRecon(x, xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint):
    

 # double x = m_mom2long/m_mom2trans;
  #//static double xEndpoint = 1.;
  #//static double yEndpoint = 1.2636363636;
  if (x < xPivot):
     yPivot = expOffset + deltaY
     m = (yPivot - yEndpoint) / (xPivot - xEndpoint)
     q = yPivot - m * xPivot
     return m * x + q
  else: 
    return deltaY * np.exp(-(x - xPivot)/expIndex) + expOffset



def ws_newRecon(x, xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint):
    

 # double x = m_mom2long/m_mom2trans;
  #//static double xEndpoint = 1.;
  #//static double yEndpoint = 1.2636363636;
  if (x < xPivot):
     yPivot = expOffset + deltaY
     m = (yPivot - yEndpoint) / (xPivot - xEndpoint)
     q = yPivot - m * xPivot
     return m * x + q
  else: 
    return deltaY * np.exp(-(x - xPivot)/expIndex) + expOffset









#def plot_Dmin(data1Dmin,data2wsDmin, data3d, base_dir):
def plot_Dmin(dataLW, base_dir):



    # PLOT: BEST Dmin vs energy

    # dati s. castellano (from scan on MC sims):
    eMC = np.array( [2.01, 2.29, 2.70, 2.98, 3.69, 4.00, 4.50, 5.00, 6.00, 6.40, 7.00, 8.00])## energy keV
    dminMC = np.array([0.3, 0.3, 0.7, 0.7, 1.7, 2.0, 2.1, 2.1, 1.9, 1.9, 1.8, 1.7])## dmin best

    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax1=plt.subplot(111)
    ax1.set_title('best dmin')

    #scan 1d
   # plt.errorbar(data1Dmin.energy,data1Dmin.best_val, fmt='bo',label='1D scan')
   # ax1.fill_between(data1Dmin.energy,data1Dmin.best_val_low, data1Dmin.best_val_up,color='gray',alpha=0.1, interpolate=True)

    
    ax1.axhline(y=1.5,label='standard_value',color='gray', linestyle='--',alpha=0.5)

    #scan 2d
    #plt.errorbar(data2wsDmin.energy,data2wsDmin.best_var1, fmt='ro',label='best dmin 2d scan')
    #ax1.fill_between(data2wsDmin.energy,data2wsDmin.best_var1_low, data2wsDmin.best_var1_up,color='red',alpha=0.1, interpolate=True)

    #scan3s=d
   # plt.errorbar(data3d.energy,data3d.best_var1, fmt='ro',label='3D scan')
   # ax1.fill_between(data3d.energy,data3d.best_var1_low, data3d.best_var1_up,color='red',alpha=0.2, interpolate=True)

    #scan in bin di LW
    plt.errorbar(dataLW.energy,dataLW.best_var1, fmt='ko',label='scan Dmin-ws')
    ax1.fill_between(dataLW.energy,dataLW.best_var1_low, dataLW.best_var1_up,color='gray',alpha=0.2, interpolate=True)



   

    # parametrizzazione
    x=np.linspace(0,10,1000)
    """
    #old:
    p0=2.4
    p1=1.15
    p2=0.03
    p3=1.87
    p4=0.63
    p5=0.5
    """
    #new:
    p0=2.05
    p1=1.15
    p2=0.03
    p3=1.87
    p4=0.3
    p5=0.645


    
    
    y=(p0*(1./(  np.exp(  -(x-p1)/p2) +1. )    ) ) * (  (1.-p5)/(np.exp( (x-p3)/p4)+1) +p5  )
    plt.plot(x,y, 'k',  label='full parametrization',alpha=0.7)

    # new parametization (piatta a bassi L/W)
    
    p0=1.8
    p1=2.2
    p2=0.11
    p3=1.33

    y2=( ((p0-p3) /(np.exp( (x-p1)/p2)+1)) +p3  )
    plt.plot(x,y2, 'g', label='new parametrization')



    #parametrizzazione nre_recon
    #--weight-scale-offset 0.04  --weight-scale-expo-idx  1.5    --weight-scale-expo-delta -0.01   --weight-scale-end-point 0.0845455  --pivot 2.5 \n

    xPivot=2.5
    deltaY=0.2
    expIndex=1.6
    expOffset=1.5 #??????????????
    xEndpoint=1
    yEndpoint =1.2636363636000001

    
    ynew=np.array([0.]*len(x))
    #print("ynew=",ynew)
    for i in range (0,len(x)):
        ynew[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)

    plt.plot(x,ynew, 'm', label='new recon')
        
   # my para 0
 # --dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.4   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0845455  --pivot 2.0 \n'                                               
   
    xPivot=2.0 
    deltaY=0.4 
    expIndex=0.6
    expOffset=1.33
    xEndpoint=1
    yEndpoint =1.2636363636000001
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
        
    #plt.plot(x,ynew_my, 'r', label='new recon-simo v0')

    ###########################################################33
    # my para 1
    
#     --dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.3   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0865455  --pivot 2.2 \n'                                          
    xPivot=2.2 
    deltaY=0.5 
    expIndex=0.6
    expOffset=1.33
    xEndpoint=1
    yEndpoint =1.2636363636000001
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
        
    #plt.plot(x,ynew_my, 'b', label='new recon-simo v1')
    #############################################################
      # my para 2
    
     #  --dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.3   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0865455  --pivot 2.5 \n'

    xPivot=2.5 
    deltaY=0.3 
    expIndex=0.6
    expOffset=1.33
    xEndpoint=1
    yEndpoint =1.2636363636000001
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
             
    #plt.plot(x,ynew_my, 'g', label='new recon-simo v2')
  
    ###################################################
    # my para 3
       
    xPivot=1.9
    deltaY=0.4 
    expIndex=0.6
    expOffset=1.33
    xEndpoint=1
    yEndpoint =1.7
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
          
    #plt.plot(x,ynew_my, 'r', label='new recon-simo v3')

 ###################################################
    # my para 5
       
    xPivot=1.9
    deltaY=0.45 
    expIndex=0.6
    expOffset=1.33
    xEndpoint=1
    yEndpoint =1.6
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
          
    plt.plot(x,ynew_my, 'c', label='new recon-simo v5')



    
    plt.xlim(1,10.5)
    plt.xlabel('elongation (M2L/M2T)')
    plt.ylabel('dmin')
    plt.legend(loc='lower right')
    #reorderLegend(ax1,['best dmin',r'1$\sigma$ band' ,'best dmin 2d scan', r'1$\sigma$ band 2d scan',  'best dmin MC', 'standard_value' ])
    #reorderLegend(ax1,['best dmin','best dmin 2d scan','3d scan' , 'best dmin MC', 'standard_value' ])
   # reorderLegend(ax1,['1D scan', '3D scan','standard_value'])
    outfilePlot=base_dir+'Dmin_scanSummary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)



###################

def  plot_Ws( dataLW, base_dir):

    # fig2# PLOT: BEST weight scale vs energy

   
    
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax2=plt.subplot(111)
    ax2.set_title('best weight scale')
    ax2.axhline(y=0.05,label='standard_value',color='gray', linestyle='--',alpha=0.5)

    #scan3d
   # plt.errorbar(data3d.energy,data3d.best_var2, fmt='ro',label='3D scan')
   # ax2.fill_between(data3d.energy,data3d.best_var2_low, data3d.best_var2_up,color='red',alpha=0.2, interpolate=True)


   #scan2d LW
    plt.errorbar(dataLW.energy,dataLW.best_var2, fmt='ko',label='scan Rmin-ws' )
    ax2.fill_between(dataLW.energy,dataLW.best_var2_low, dataLW.best_var2_up,color='gray',alpha=0.2, interpolate=True)

    #plt.plot(eMC_ws,wsMC,'ko',mfc='none',markersize=10, label='best ws MC')
   # parametrizzazione
    x=np.linspace(0,10,1000)
    """
    p0=0.341213
    p1=1.21184
    p2=0.0824731
    p3=0.03
    y=p0*np.exp(-0.5*( ((x-p1)/p2) + np.exp( (p1-x)/p2)    )   )+p3
    """
        
    p0= 0.341213 
    p1= 1.21184
    p2= 0.0824731
    p3= 0.028 
    p4= 0.008
    p5= 5.5
    p6= 1

    y=p0*np.exp(-0.5*( ((x-p1)/p2) + np.exp( (p1-x)/p2)    )   )+p3+ (  ( ((p4) /(np.exp( -(x-p5)/p6)+1))  )    ) 
    
    plt.plot(x,y,'k', label='parametrization full')


    # new parametization (piatta a bassi L/W)
    
    p0=0.1
    p1=1.57
    p2=0.1
    p3=0.03

    y2=( ((p0-p3) /(np.exp( (x-p1)/p2)+1)) +p3  )
   # plt.plot(x,y2,'b', label='new parametrization (1)')

    p3=0.028
    p4= 0.006 
    p5= 4.7
    p6= 0.512
  
    y3=( ((p0-p3) /(np.exp( (x-p1)/p2)+1)) +p3  ) +  ( ((p4) /(np.exp( -(x-p5)/p6)+1))  )
    plt.plot(x,y3,'g', label='new parametrization (2)')


    

    #parametrizzazione ixperecon nuova s.c.

    #parametrizzazione nre_recon
    #--weight-scale-offset 0.04  --weight-scale-expo-idx  1.5    --weight-scale-expo-delta -0.01   --weight-scale-end-point 0.0845455  --pivot 2.5 \n

   # parametri sc
    xPivot=2.5
    deltaY=-0.01
    expIndex=1.5
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.0845455

    ynew=np.array([0.]*len(x))
    #print("ynew=",ynew)
    for i in range (0,len(x)):
        ynew[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
        
    plt.plot(x,ynew, 'm', label='new recon')
    
    ################################################
     # my para 0
     # --dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.4   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0845455  --pivot 2.0 \n'                                               
   
    xPivot=2.0
    deltaY=-0.012
    expIndex=10
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.0845455
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
 #   plt.plot(x,ynew_my, 'r', label='new recon - my mod v0')
    
    ###########################################################33
    # my para 1
    
    #--dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.3   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0865455  --pivot 2.2 \n'   


    xPivot=2.2
    deltaY=-0.012
    expIndex=10
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.0865455
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
  #  plt.plot(x,ynew_my, 'b', label='new recon - my mod v1')

    ##################################################
    # my para 2
    
    #  --dmin-offset 1.33 --dmin-expo-idx 0.6    --dmin-expo-delta 0.3   --dmin-end-point 1.2636363636000001      --weight-scale-offset 0.04  --weight-scale-expo-idx  10.    --weight-scale-expo-delta -0.012 --weight-scale-end-point 0.0865455  --pivot 2.5 \n'

    xPivot=2.5
    deltaY=-0.012
    expIndex=10
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.0865455
    
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
  #  plt.plot(x,ynew_my, 'g', label='new recon - my mod v2')
     


     ##################################################
    # my para 3
    
    xPivot=1.9
    deltaY=-0.012
    expIndex=10
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.2

        
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
   # plt.plot(x,ynew_my, 'g', label='new recon - my mod v3')
     
 ##################################################
    # my para 5
    
    xPivot=1.9
    deltaY=-0.011
    expIndex=10
    expOffset=0.04
    xEndpoint=1
    yEndpoint =0.22

        
    ynew_my=np.array([0.]*len(x))
    for i in range (0,len(x)):
        ynew_my[i]=dmin_newRecon(x[i], xPivot, deltaY,  expIndex,  expOffset,  xEndpoint,  yEndpoint)
        
    plt.plot(x,ynew_my, 'c', label='new recon - my mod v5')
     



    ##########################

    
   # plt.xlim(1,7)
    plt.xlabel('elongation (M2L/M2T)')
    plt.ylabel('weight scale')
    plt.legend()
    
    #reorderLegend(ax2,['best ws','best ws 2d scan', '3d scan',  'best ws MC', 'standard_value' ])
    reorderLegend(ax2,['1D scan', '3D scan','standard_value'])
    outfilePlot=base_dir+'ws_scanSummary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)




def plot_Dmax(data1Dmin,data2wsDmin, data3d, base_dir):
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
    
    #plt.plot(eMC_dmax,dmaxMC,'ko',mfc='none',markersize=10, label='best Dmax MC')
    plt.xlabel('energy [KeV]')
    plt.ylabel('best Dmax')
    plt.legend(loc='upper left')
    #reorderLegend(ax4,['best Dmax', 'best Dmax MC', 'standard_value' ])



    outfilePlot=base_dir+'dmax_scanSummary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)

#######################



##############################################
def plot_Moma1(data1Dmin,data2wsDmin, data3d, base_dir):

    # PLOT moma1

    # dati scan MC(s.castellano)
    moma1MC =np.array( [22,28,34,32,22,20,20])#moma1 best (energie tue senza aggiunte)

   
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.9, bottom=0.09,hspace=0.250)
    ax5=plt.subplot(111)
    ax5.set_title('best moma1')
    plt.errorbar(data1Moma1.energy,data1Moma1.best_val, fmt='bo',label='1D scan')
    ax5.fill_between(data1Moma1.energy,data1Moma1.best_val_low, data1Moma1.best_val_up,color='gray',alpha=0.1, interpolate=True)
    ax5.axhline(y=36,label='standard_value', linestyle='--',alpha=0.5)


    #plt.errorbar(data2zeroMoma.energy,data2zeroMoma.best_var2, fmt='ro',label='best moma1(2) 2D scan')
    #ax5.fill_between(data2zeroMoma.energy,data2zeroMoma.best_var2_low, data2zeroMoma.best_var2_up,color='red',alpha=0.1, interpolate=True )

    #plt.plot(data1Moma1.energy,moma1MC,'ko',mfc='none',markersize=10, label='best moma1 MC')


    #scan 3d
    plt.errorbar(data3d.energy,data3d.best_var3, fmt='ro',label='3D scan')
    ax5.fill_between(data3d.energy,data3d.best_var3_low, data3d.best_var3_up,color='red',alpha=0.2, interpolate=True)


    # parametrizzazione
    x=np.linspace(0,8,100)
    p0= 36
    p1=20
    p2=3.33
    p3=0.08

    y3=  (p0-p1)/(np.exp( (x-p2)/p3)+1) +p1 

 #   plt.plot(x,y3, label='parametrization')

    
    plt.xlim(1,7)
    plt.xlabel('energy [KeV]')
    plt.ylabel('best moma1')
    plt.legend()
    reorderLegend(ax5,['1D scan', '3D scan','standard_value'])


    outfilePlot=base_dir+'moma1_scanSummary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)



#########################################

def plot_Moma2(data1Dmin,data2wsDmin, data3d, base_dir):
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
    #reorderLegend(ax6,['best Moma2', 'best moma1(2) 2D scan','best moma2 MC', 'standard_value' ])

    outfilePlot=base_dir+'moma2_scanSummary.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)


#########################################
def plot_ClusteringThr(data1Dmin,data2wsDmin, data3d, base_dir):

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

def plot_Modulation(data1Dmin,data2wsDmin, data3d, dataBest, base_dir):

    # PLOT modulation factor

   
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.093, right=0.97, top=0.90, bottom=0.09,hspace=0.34,wspace=0.255)
    ax3=plt.subplot(111)
    
    plt.grid(True,linestyle=':', color='grey') 

    ax3.set_title('modPhi2best/ max(modPhi1std,modPhi2std)')
    plt.errorbar(data1ws.energy,data1ws.mod2/(np.maximum(data1ws.mod2std,data1ws.mod1std)), fmt='bo-',label='scan ws')
    plt.errorbar(data1Dmin.energy,data1Dmin.mod2/(np.maximum(data1Dmin.mod2std,data1Dmin.mod1std)), fmt='ro-',label='scan dmin')
  #  plt.errorbar(data2wsDmin.energy,data2wsDmin.mod2/(np.maximum(data2wsDmin.mod2std,data2wsDmin.mod1std)), fmt='ko',label='scan dmin-ws')

    #plt.errorbar(data1ZeroThr.energy, (np.maximum(data1ZeroThr.mod2, data1ZeroThr.mod1 ))/(np.maximum(data1ZeroThr.mod2std,data1ZeroThr.mod1std)), fmt='go-',label='scan clustering thr (max mo1d mod2)')
    plt.errorbar(data1Moma1.energy,(np.maximum(data1Moma1.mod2,data1Moma1.mod1  ))/(np.maximum(data1Moma1.mod2std,data1Moma1.mod1std)), fmt='mo-',label='scan moma1 thr')
    
    plt.errorbar(data1Moma2.energy,data1Moma2.mod2/(np.maximum(data1Moma2.mod2std,data1Moma2.mod1std)), fmt='ko-',label='scan moma2 thr')
    plt.errorbar(data1Dmax.energy,data1Dmax.mod2/(np.maximum(data1Dmax.mod2std,data1Dmax.mod1std)), fmt='go-',label='scan dmax')

    #plt.errorbar(data2zeroMoma.energy,(np.maximum(data2zeroMoma.mod2,data2zeroMoma.mod1 ))/(np.maximum(data2zeroMoma.mod2std,data2zeroMoma.mod1std)), fmt='kP',label='scan zero-moma12  (max mod1 mod2)  ')


    #plt.errorbar(data3d.energy,data3d.mod2/np.maximum(data3d.mod2std,data3d.mod1std), fmt='mo',label='scan 3D Dmin-Ws - moma1,2  ')

  #  plt.errorbar(dataBest.energy,dataBest.mod2/np.maximum(dataBest.mod2std,dataBest.mod1std), fmt='mP--',label='scan Dmin-Ws - moma1,2 - NEW DATA  ')
    #plt.errorbar(dataBest.energy,dataBest.mod1/np.maximum(dataBest.mod2std,dataBest.mod1std), fmt='m*--',label='scan Dmin-Ws - moma1,2 - PHI1   NEW DATA  ')

    plt.xlabel('energy [KeV]')
    plt.ylabel('modPhi2best/ max(modPhi1std,modPhi2std)')
    plt.legend()
    outfilePlot=base_dir+'modComparison.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)

###################################33


def plot_ModulationAll( dataBest_reconPara,  dataBest_reconParaLW_new1,  dataBest_reconParaLW_new2, dataBest_reconParaLW_old, dataBest_ixperecon_new,dataBest_ixperecon_my, dataBest_ixperecon_myV2, dataBest_ixperecon_myV3,  dataBest_ixperecon_myV4, dataBest_ixperecon_myV5,dataBest_reconParaLW_FullUpdate, base_dir):

     
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.093, right=0.97, top=0.90, bottom=0.09,hspace=0.34,wspace=0.255)
    ax=plt.subplot(111)
    plt.grid(True,linestyle=':', color='grey') 
    ax.set_title('mod Phi2 best/mod. std ( max phi1, phi2)')

    #mod_stdAll=unisci_phi(dataBest.energy,dataBest_reconPara.mod1std,dataBest.mod2std)

    #mod_stdAll_scan=unisci_phi(data3d.energy,data3d.mod1std,data3d.mod2std)
 


    modStd_max=np.maximum(dataBest_reconPara.mod2std,dataBest_reconPara.mod1std)
    
    
    plt.errorbar( dataBest_reconPara.energy, dataBest_reconPara.mod2/modStd_max  , fmt='ko--',label='parametrizzazione PHA', alpha=0.4)
#    plt.errorbar( dataBest_reconParaLW_new1.energy, dataBest_reconParaLW_new1.mod2/modStd_max  , fmt='bo-',label='parametrizzazione M2L/M2T new1', alpha=1)
    plt.errorbar( dataBest_reconParaLW_new2.energy, dataBest_reconParaLW_new2.mod2/modStd_max, fmt='go-',label='parametrizzazione M2L/M2T new2', alpha=1)
   # plt.errorbar( dataBest_reconParaLW_old.energy, dataBest_reconParaLW_old.mod2/modStd_max , fmt='ko-',label='parametrizzazione M2L/M2T old', alpha=1)
   
   

    plt.errorbar( dataBest_ixperecon_new.energy, dataBest_ixperecon_new.mod2/modStd_max , fmt='ms-',label='parametrizzazione branch ixperecon_new', alpha=1)
   
    #plt.errorbar( dataBest_ixperecon_my.energy, dataBest_ixperecon_my.mod2/modStd_max , fmt='rs-',label='parametrizzazione branch ixperecon_new- modificata v0', alpha=1)
   
  #  plt.errorbar( dataBest_ixperecon_myV2.energy, dataBest_ixperecon_myV2.mod2/modStd_max , fmt='bP-',label='parametrizzazione branch ixperecon_new- modificata V1', alpha=1)
   
  #  plt.errorbar( dataBest_ixperecon_myV3.energy, dataBest_ixperecon_myV3.mod2/modStd_max , fmt='gX-',label='parametrizzazione branch ixperecon_new- modificata V2', alpha=1)

   # plt.errorbar( dataBest_ixperecon_myV4.energy, dataBest_ixperecon_myV4.mod2/modStd_max , fmt='yX-',label='parametrizzazione branch ixperecon_new- modificata V4', alpha=1)

    plt.errorbar( dataBest_ixperecon_myV5.energy, dataBest_ixperecon_myV5.mod2/modStd_max , fmt='cX-',label='parametrizzazione branch ixperecon_new- modificata V5', alpha=1)
 
    plt.errorbar( dataBest_reconParaLW_FullUpdate.energy, dataBest_reconParaLW_FullUpdate.mod2/modStd_max , fmt='kX-',label='parametrizzazione Full updated', alpha=1)
 
    plt.xlabel('energy [KeV]')
    #plt.ylabel('mod Phi2_best/ mod. std [max(phi1,phi2)]')

    plt.legend()
    outfilePlot=base_dir+'modComparisonAll_LW.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)               



    

    
    
########################
def plot_Phase(data1Dmin,data2wsDmin, data3d, dataBest,dataBest_reconPara ,  base_dir):
   

   # fare con dataset indipendente!!!!!!!!!!!!!!!!!!!!!! TO DO
    
    # plot phase best-std
    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.093, right=0.97, top=0.90, bottom=0.09,hspace=0.34,wspace=0.255)
    ax8=plt.subplot(111)
    plt.grid(True,linestyle=':', color='grey') 

    ax8.set_title('phase')


    delta2=data3d.phase2-data3d.phase2std
    delta2_err=(data3d.phase2_err**2+data3d.phase2std_err**2)**0.5 
    delta1=data3d.phase1-data3d.phase1std
    delta1_err=(data3d.phase1_err**2+data3d.phase1std_err**2)**0.5 


    #plt.errorbar(data3d.energy,delta2,yerr=delta2_err, fmt='bo',label='delta phase2')
    #plt.errorbar(data3d.energy, delta1,yerr=delta1_err, fmt='ro',label='delta phase1')

    phase_std=unisci_phi(dataBest_reconPara.energy,dataBest_reconPara.phase1std,dataBest_reconPara.phase2std)
    phase_std_err=unisci_phi(dataBest_reconPara.energy,dataBest_reconPara.phase1std_err,dataBest_reconPara.phase2std_err)
    
    
    plt.errorbar(dataBest_reconPara.energy,  phase_std  ,yerr=phase_std_err, fmt='ro',label='std')

    plt.errorbar(dataBest_reconPara.energy,dataBest_reconPara.phase2 ,yerr=dataBest_reconPara.phase2_err, fmt='bo',label='PHA parametrization')


    plt.xlabel('energy [KeV]')
    plt.ylabel('phase [deg]')
    plt.legend()

    outfilePlot=base_dir+'Phase.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)


    fig=plt.figure(figsize=(10,7))
    fig.subplots_adjust(left=0.093, right=0.97, top=0.90, bottom=0.09,hspace=0.34,wspace=0.255)
    ax9=plt.subplot(111)
    plt.grid(True,linestyle=':', color='grey') 
    ax9.set_title('phase best - phase std')

    plt.errorbar(dataBest_reconPara.energy, dataBest_reconPara.phase2- phase_std,   fmt='ro')


    
    outfilePlot=base_dir+'deltaPhase.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)


#def creaRootFiles(data1Dmin,data2wsDmin, data3d,dataBest, base_dir):
def creaRootFiles(dataLW, base_dir):

    
    # crea root TGraph per paratmetri ottimali
    import ROOT
    from array import array


    gDmin_lw=ROOT.TGraph(len(dataLW.energy),array('f', dataLW.energy.tolist()),array('f', dataLW.best_var1.tolist()) )
    gDmin_lw.SetName('gDmin_lw')
    gDmin_lw.SetTitle('best Dmin_lw vs L/W')
    gDmin_lw.SetMarkerStyle(20)
    gDmin_lw.SetMarkerColor(2)
    gDmin_lw.Draw("ap")

    errUp=dataLW.best_var1_up-dataLW.best_var1
    errLow=dataLW.best_var1-dataLW.best_var1_low
    errX=array('f',[0.]*len(dataLW.energy))

    gDmin_lwAsymErr=ROOT.TGraphAsymmErrors(len(dataLW.energy),array('f', dataLW.energy.tolist()),array('f', dataLW.best_var1.tolist()),errX,errX,array('f', errLow.tolist()), array('f', errUp.tolist()))
    gDmin_lwAsymErr.SetName('gDmin_lwAsymm')
    gDmin_lwAsymErr.Draw("*")
   

    
    # plot ws vs LW
    gWs_lw=ROOT.TGraph(len(dataLW.energy),array('f', dataLW.energy.tolist()),array('f', dataLW.best_var2.tolist()) )
    gWs_lw.SetName('gWs_lw')
    gWs_lw.SetTitle('best Ws vs LW')
    gWs_lw.SetMarkerStyle(20)
    gWs_lw.SetMarkerColor(2)
    gWs_lw.Draw("ap")

    
       
    outRootFile=ROOT.TFile('/home/maldera/IXPE/rec_optimization/parametrizzazione_LW.root','recreate')
    gDmin_lw.Write()
    gDmin_lwAsymErr.Write()
    gWs_lw.Write()

    #gMoma_3d.Write()
    #gDmin_3dPHA.Write()
    #gWs_3dPHA.Write()
    #gMoma_3dPHA.Write()
    outRootFile.Close()




def plotVsPha(data1Dmin,data2wsDmin, data3d, dataBest,  base_dir):
    mpl.rcParams['font.size']=14  #!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    fig=plt.figure(figsize=(10,8))
    fig.subplots_adjust(left=0.085, right=0.97, top=0.95, bottom=0.09,hspace=0.34,wspace=0.255)

    x=dataBest.pha
    xPara=np.linspace(6000,20000,100)

    plt.grid(True,linestyle=':', color='grey') 
    #moma12 
    ax=plt.subplot(221)
    ax.set_title('moma1,2 threshold')
    plt.errorbar(x,data3d.best_var3, fmt='ro',label='3d scan')
    ax.fill_between(x,data3d.best_var3_low, data3d.best_var3_up,color='red',alpha=0.1, interpolate=True)
    ax.axhline(y=36,label='standard_value', linestyle='--',alpha=0.5)
    #parametrizzazione:
    p0= 36
    p1=20
    p2=10339.5
    p3=250

    y3=  (p0-p1)/(np.exp( (xPara-p2)/p3)+1) +p1 

    plt.plot(xPara,y3, label='parametrization')
    plt.xlabel('pha',x=0.9)
    plt.ylabel('moma threshold')
    plt.legend(loc='upper right')




    #Dmin 
    ax=plt.subplot(222)
    ax.set_title('Dmin')
    
    ax.axhline(y=1.5,label='standard_value', linestyle='--',alpha=0.5)
    plt.errorbar(x,data3d.best_var1, fmt='ro',label='3d scan')
    ax.fill_between(x,data3d.best_var1_low, data3d.best_var1_up,color='red',alpha=0.1, interpolate=True)

    p0=1.5406
    p1=8107.25
    p2=250
    p3=11972.5
    p4=770.5
    p5=0.787144

    y=(p0*(1./(  np.exp(  -(xPara-p1)/p2) +1. )    )   +0.2) * (  (1.-p5)/(np.exp( (xPara-p3)/p4)+1) +p5  )


    plt.plot(xPara,y, label='parametrization')
   
    plt.xlabel('pha',x=0.9)
    plt.ylabel('Dmin')
    plt.legend()

 

    
    #Ws 
    ax=plt.subplot(223)
    ax.set_title('Weight scale')
    plt.errorbar(x,data3d.best_var2, fmt='ro',label='3d scan')
    ax.fill_between(x,data3d.best_var2_low, data3d.best_var2_up,color='red',alpha=0.1, interpolate=True)

    ax.axhline(y=0.05,label='standard_value', linestyle='--',alpha=0.5)

    p0=0.352
    p1=0.03
    p2=8551
    p3=813
    p4=0.31
    p5=7200
    p6=380
    
    y= ((p0-p1)/(np.exp( (xPara-p2)/p3)+1) +p1  ) *  (( ( 1-p4) / ( np.exp( -(xPara-p5)/p6) +1. ))+p4)

    plt.plot(xPara,y, label='parametrization')
    plt.xlabel('pha',x=0.9)
    plt.ylabel('ws')
    plt.legend()

    
    
    
    #ratio vs E
   # ax=plt.subplot(224)
    #ax.set_title('PHA/Energy vs E ')
    #plt.errorbar(data3d.energy, dataBest.pha/dataBest.energy , fmt='go',label='')
   

    outfilePlot=base_dir+'pha_parametrization.png'
    print ("outFile png =",outfilePlot)
    plt.savefig(outfilePlot)
    
    
    

######################

base_dir='/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization//scanRmin-Ws/LWbins/'



data1Dmin=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin/outScan_dmin_v2.txt')
data1ws=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanWs/outScan_weight_scale_v2.txt')

data2wsDmin=readData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/outScan_dmin-weight_scale.txt')


#tutti i plot per rapporto 
data1ZeroThr=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/outScan_zero_thr_v2.txt')
data1Moma1=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/outScan_moma1_thr_v2.txt')
data1Moma2=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/outScan_moma2_thr_v2.txt')
data1Dmax=readSummaryData_1d_v2('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmax/outScan_dmax_v2.txt')
data2zeroMoma=readData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanZeroMoma12/outScan_zero_thr-moma1_thr.txt')


# data 3D

data3d=readData3d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/outScan_dmin-weight_scale_v2.txt')





# data set indipendente!!
dataBest=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_v2.txt', version=2)

dataBest_para=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_para_v2.txt', version=2)
dataBest_reconPara=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixperecon_para_v2.txt', version=2)


dataBest_reconParaLW_new1=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixpereconPara_LW_new1_v2.txt', version=2)
dataBest_reconParaLW_new2=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixpereconPara_LW_new2_v2.txt', version=2)
dataBest_reconParaLW_old=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixpereconPara_LW_v2.txt', version=2)

dataBest_ixperecon_new=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixprecon_new_v2.txt', version=2)

dataBest_ixperecon_my=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_my_ixperecon_para_v2.txt', version=2)
dataBest_ixperecon_myV2=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_my_ixperecon_para_v2_v2.txt', version=2)
dataBest_ixperecon_myV3=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_my_ixperecon_para_v3_v2.txt', version=2)


dataBest_ixperecon_myV4=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_my_ixperecon_para_v4_v2.txt', version=2)

dataBest_ixperecon_myV5=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_my_ixperecon_para_v5_v2.txt', version=2)

dataBest_reconParaLW_FullUpdate=readDataBest('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/bestParams/outBestParams_ixpereconPara_LW_fullUpdated_v2.txt', version=2)



# scan in LW

dataLW=readData2d('/home/maldera/IXPE/rec_optimization/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/LWbins/outScanLW_dmin-weight_scale.txt')



######################
# print obtimal values:

print("\n \n BEST PARAMETERS")
print("bestZeroThr=",data2zeroMoma.best_var1)
print("bestMoma12=",data2zeroMoma.best_var2)
print("bestDmin=",data2wsDmin.best_var1)
print("bestws=",data2wsDmin.best_var2)
print("E=",data2wsDmin.energy)


print("\n \n BEST PARAMETERS  3D  ")
#print("bestZeroThr=",data2zeroMoma.best_var1)
print("bestMoma12=",data3d.best_var3)
print("bestDmin=",data3d.best_var1)
print("bestws=",data3d.best_var2)
print("E=",data3d.energy)



#plot_Dmin(data1Dmin,data2wsDmin, data3d, base_dir)
plot_Dmin(dataLW, base_dir)




#plot_Ws(data1Dmin,data2wsDmin, data3d, base_dir)
plot_Ws(dataLW, base_dir)



#plot_Dmax(data1Dmin,data2wsDmin, data3d, base_dir)

#plot_Moma1(data1Dmin,data2wsDmin, data3d, base_dir)

#plot_Moma2(data1Dmin,data2wsDmin, data3d, base_dir)


#plot_ClusteringThr(data1Dmin,data2wsDmin, data3d, base_dir)



#plot_Modulation(data1Dmin,data2wsDmin, data3d, dataBest, base_dir)


#plot_Phase(data1Dmin,data2wsDmin, data3d, dataBest,dataBest_reconPara, base_dir)

creaRootFiles(dataLW, base_dir)



#plotVsPha(data1Dmin,data2wsDmin, data3d, dataBest, base_dir)

#plot_ModulationAll(data1Dmin,data2wsDmin, data3d, dataBest, dataBest_para, dataBest_reconPara,   base_dir)
plot_ModulationAll(dataBest_reconPara,  dataBest_reconParaLW_new1,  dataBest_reconParaLW_new2, dataBest_reconParaLW_old, dataBest_ixperecon_new , dataBest_ixperecon_my,dataBest_ixperecon_myV2,  dataBest_ixperecon_myV3,  dataBest_ixperecon_myV4  ,  dataBest_ixperecon_myV5, dataBest_reconParaLW_FullUpdate,  base_dir)


plt.show()
