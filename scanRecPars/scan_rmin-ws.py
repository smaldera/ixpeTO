import lanciaRecon_hct
import glob
import numpy as np


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


base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/'
dirs=['001361',  '001388',  '001416',  '001436',  '001461',  '001471','001333']
#dirs=['001361',  '001388']


for dir in dirs:
    dataDir='/data1/IXPE/data/misureDU_2/'+dir
    print("ecco i files:")
    n_files=1
    if dir=='001461' or dir=='001471':
       n_files=10

    filenames=get_files(dataDir,n_files)
    print ("primi 4 files= ",filenames)
 
    out_dir=base_dir+dir+'/'


    n_iter=187
    #inizio scan su zero_thr
    #for dmin  in np.arange(1.1,1.5,0.2):
    for dmin  in [1.6, 1.8,2.0]: # scan fine     
      for weight_scale  in np.arange(0.03,0.25,0.02):
        #for weight_scale  in [0.029,0.031]:
        
       
          moma2_thr=moma1_thr        
          print (n_iter," dmin = ",dmin, " weight_scale =",weight_scale)
          work_dir=out_dir+str(n_iter)
          
          lanciaRecon_hct.submit_recon(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir)

          filename=work_dir+'/config_simo.txt'
          myLogfile=open(filename,'w')  
          myLogfile.write("zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
          
          n_iter=n_iter+1
          


    
    
    

