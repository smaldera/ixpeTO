import lanciaRecon_hct
import glob



def get_files(dataDir,n,n_max):

 nomi_out=''
 for i in range (1,n_max+1):
    nomi_out+=' '+dataDir+str(i)+'/out_sim.fits'
    if i>n: break

 
 #f=glob.glob(dataDir+"/*_1118.lv1.fits")
 ##print(f)
 #creo stringa con i primi n files
 #i=0
 #nomi_out=''
# for name in f:
#     nomi_out+=" "+name
#     if i>(n-1): break
#     i+=1
 return nomi_out




#filenames="/data1/IXPE/data/misureDU_2/001333/ixpe2019100*_1118.lv1.fits"
filenames="/data1/IXPE/data/misureDU_2/001333/*.1??_1203_001333_1118.lv1.fits"  # 4 FILES!!!
output_folder='.'
zero_thr=15
moma1_thr=36
moma2_thr=36
coer_noise_offset=0
trig_minicluster_offset=0
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


out_base_dir='/data1/maldera/IXPE_work/rec_optimization/simData/scan_mona2Thr/'
data_base_dir='/data1/maldera/IXPE_work/rec_optimization/simData/'

#dirs=['2KeV',  '3KeV',  '4KeV',  '5KeV',  '6KeV',  '7KeV','8KeV']
#dirs=['2KeV',  '3KeV',  '4KeV',   '6KeV',  '7KeV','8KeV']
dirs=[ '5KeV']


angle_dirs=['angle-180deg']
max_dirs={'2KeV':29,  '3KeV':17,  '4KeV':0,  '5KeV':5,  '6KeV':3,  '7KeV':29,'8KeV':18}


for dir in dirs:
  for  angle_dir in angle_dirs:
    dataDir=data_base_dir+dir+'/'+angle_dir+'/'
    
    print("ecco i files:")
    n_files=max_dirs[dir]
    
    filenames=get_files(dataDir,n_files,max_dirs[dir])
    print ("files= ",filenames)
 
    out_dir=out_base_dir+dir+'/'+angle_dir+'/'
    n_iter=1
    #inizio scan su zero_thr
    for moma2_thr in range (16,45,2):

        print ("moma2_thr= ",moma2_thr)
        work_dir=out_dir+str(n_iter)
        
        lanciaRecon_hct.submit_recon(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir)

        filename=work_dir+'/config_simo.txt'
        myLogfile=open(filename,'w')  
        myLogfile.write("zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
        n_iter=n_iter+1



    
    
    

