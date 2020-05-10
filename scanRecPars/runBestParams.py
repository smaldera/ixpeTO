import lanciaRecon_hct
import glob

import numpy as np


"""
def get_files(dataDir,n):

 f=glob.glob(dataDir+"/*_1118.lv1.fits")
 #print(f)
 #creo stringa con i primi n files
 i=0
 nomi_out=''
 for name in f:
     nomi_out+=" "+name
     if i>(n-1): break
     i+=1
 return nomi_out
"""


def get_files(dataDir,n, first):
      # return n files starting fom first
      # utile per cambiare dataset velocemente...

      f=glob.glob(dataDir+"/*_1118.lv1.fits")
      #print(f)
      #creo stringa con i primi n files

      nomi_out=''
      
      for i in range(first,first+n):
            nomi_out+="  "+f[i]
            print("i=",i," f = ",f[i])
            
      return nomi_out
                                       


def param_Dmin(E):
      p0=1.58
      p1=2.67
      p2=0.10
      p3=3.753
      p4=0.2788
      p5=0.7684

      dmin=(p0*(1./(  np.exp(  -(E-p1)/p2) +1. )    )   +0.2) * (  (1.-p5)/(np.exp( (E-p3)/p4)+1) +p5)  

      return dmin                                                             
                                                                   

                                                                   

def param_Ws(E):
      p0= 0.859
      p1=0.03
      p2=2.658
      p3=0.21
      p4=0.1
      p5=2.7
      p6=0.246
      
      ws=( p0*(1.-p1)/(np.exp( (E-p2)/p3)+1) +p1  )  *  (( ( 1-p4) / ( np.exp( -(E-p5)/p6) +1. ))+p4)
      return ws

def param_Moma12(E):
    p0= 36
    p1=20
    p2=3.33
    p3=0.08

    y3=  (p0-p1)/(np.exp( (E-p2)/p3)+1) +p1 

    moma12=round(y3)
    return moma12




#filenames="/data1/IXPE/data/misureDU_2/001333/ixpe2019100*_1118.lv1.fits"
filenames="/data1/IXPE/data/misureDU_2/001333/*.1??_1203_001333_1118.lv1.fits"  # 4 FILES!!!
output_folder='.'
zero_thr=20
moma1_thr=36
moma2_thr=36
coer_noise_offset=3
trig_minicluster_offset=2
suffix='recon'
min_track_hits=6
min_densy_points=4
dmin=1.5
dmax=3.5
weight_scale=0.05
truncate_lsb=0
logfile="recon.log"
n_max=1000000000
first_ev=0


#E= [2.01 2.29 2.7  2.98 3.69 4.5  6.4 ]
#bestZeroThr= {2.01:22., 2.29:20.,   2.7:20.,  2.98:22.,  3.69:20.,   4.5:18.,   6.4:12.}
bestMoma12=  {2.01:36., 2.29:36.,   2.7:36.,  2.98:36.,  3.69:20.,   4.5:20.,   6.4:20.}
bestDmin=    {2.01:0.3, 2.29:0.1,   2.7:1.1,  2.98:1.7,  3.69:1.6,   4.5: 1.4,  6.4:1.3}
bestws=      {2.01:0.11,2.29:0.19,  2.7:0.23, 2.98:0.13, 3.69:0.04,  4.5: 0.03, 6.4:0.03}



stdMoma12={2.01:36., 2.29:36.,   2.7:36.,  2.98:36.,  3.69:36.,   4.5:36.,   6.4:36.}
stdDmin={2.01:1.5, 2.29:1.5,   2.7:1.5,  2.98:1.5,  3.69:1.5,   4.5: 1.5,  6.4:1.5}
stdws={2.01:0.05,2.29:0.05,  2.7:0.05, 2.98:0.05, 3.69:0.05,  4.5: 0.05, 6.4:0.05}



paraMoma12={2.01:param_Moma12(2.01), 2.29:param_Moma12(2.29),   2.7:param_Moma12(2.7),  2.98:param_Moma12(2.98),  3.69:param_Moma12(3.69),   4.5:param_Moma12(4.5),   6.4:param_Moma12(6.4)}

paraDmin={2.01:param_Dmin(2.01), 2.29:param_Dmin(2.29),   2.7:param_Dmin(2.7),  2.98:param_Dmin(2.98) ,  3.69:param_Dmin(3.69) ,   4.5:param_Dmin(4.5) ,  6.4:param_Dmin(6.4)}


paraWs={2.01:param_Ws(2.01), 2.29:param_Ws(2.29),   2.7:param_Ws(2.7),  2.98:param_Ws(2.98),  3.69:param_Ws(3.69),   4.5:param_Ws(4.5),  6.4:param_Ws(6.4)}






dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV

base_dir='/data1/maldera/IXPE_work/rec_optimization/bestParams/'
dirs=['001361',  '001388',  '001416',  '001436',  '001461',  '001471','001333']

#dirs=['001461',  '001471','001333']


#parameters_sets=['best3d','std','para']
parameters_sets=['ixperecon_para']



dict_moma12={'best3d':bestMoma12, 'std':stdMoma12, 'para':paraMoma12}
dict_dmin={'best3d':bestDmin, 'std':stdDmin, 'para':paraDmin}
dict_ws={'best3d':bestws, 'std':stdws, 'para':paraWs}



print("parametri PARA, Moma",dict_moma12['para'])
print("parametri PARA, Ws",dict_ws['para'])
print("parametri PARA, dmin",dict_dmin['para'])





for dir in dirs:

     
    dataDir='/data1/IXPE/data/misureDU_2/'+dir
    print("ecco i files:")
    first_file=4+1  # per lo scan era implicitamente 0
    n_files=4+1
    if dir=='001461':# or dir=='001471':
       n_files=5   # 16 files in tutto!
       first_file=10+1 # nello scan era implicitaente 0

    if dir=='001471':
          n_files=4 # 15 files in tutto
          first_file=10+1 # nello scan era implicitaente 0

                             
       
    filenames=get_files(dataDir,n_files,first_file)
    print ("primi 4 files= ",filenames)


    out_dir=base_dir+dir+'/'


    for parameters in parameters_sets:

          
 
          energy=dict_energy[dir]

          if parameters != 'ixperecon_para':  #in questo caso la parametrizzaione e' gia' in ixperecon, e i dictionary non hanno questa chiave
                moma1_thr=dict_moma12[parameters][energy]
                moma2_thr=dict_moma12[parameters][energy]
                dmin=dict_dmin[parameters][energy]
                weight_scale=dict_ws[parameters][energy]
                
          else:
                print("USING std parametes, ma con parametrizzazione in ixperecon ")
                
          print("dir=",dir," E = ",dict_energy[dir]," zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
          work_dir=out_dir+parameters
        
          lanciaRecon_hct.submit_recon_local(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir)



    
          filename=work_dir+'/config_simo.txt'
          myLogfile=open(filename,'w')  
          myLogfile.write("zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )



    
    
    

