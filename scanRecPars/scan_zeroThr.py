import lanciaRecon_hct



#filenames="/data1/IXPE/data/misureDU_2/001333/ixpe2019100*_1118.lv1.fits"
filenames="/data1/IXPE/data/misureDU_2/001333/*.1??_1203_001333_1118.lv1.fits"  # 4 FILES!!!
output_folder='.'
zero_thr=30
moma1_thr=26
moma2_thr=26
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
n_max=1e10
first_ev=0


base_dir='/data1/maldera/IXPE_work/rec_optimization/scanZeroThr/'

n_iter=1
#inizio scan su zero_thr
for zero_thr in range (5,40,2):

    print ("zeroThr= ",zero_thr)
    work_dir=base_dir+str(n_iter)
    lanciaRecon_hct.submit_recon(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir)

    filename=work_dir+'/config_simo.txt'
    myLogfile=open(filename,'w')  
    myLogfile.write("zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
    n_iter=n_iter+1



    
    
    

