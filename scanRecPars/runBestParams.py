import lanciaRecon_hct
import glob



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


#E= [2.01 2.29 2.7  2.98 3.69 4.5  6.4 ]
bestZeroThr= {2.01:22., 2.29:20.,   2.7:20.,  2.98:22.,  3.69:20.,   4.5:18.,   6.4:12.}
bestMoma12=  {2.01:28., 2.29:32.,   2.7:26.,  2.98:22.,  3.69:20.,   4.5:18.,   6.4:12.}
bestDmin=    {2.01:0.3, 2.29:0.1,   2.7:1.1,  2.98:1.7,  3.69:1.9,   4.5: 1.5,  6.4:1.3}
bestws=      {2.01:0.11,2.29:0.19,  2.7:0.23, 2.98:0.13, 3.69:0.05,  4.5: 0.03, 6.4:0.03}

dict_energy={'001333':6.40, '001361':4.50,  '001388':2.98,  '001416':2.70,  '001436':2.29,  '001461':2.01,  '001471':3.69} # Energy in KeV



base_dir='/data1/maldera/IXPE_work/rec_optimization/bestParams/'
dirs=['001361',  '001388',  '001416',  '001436',  '001461',  '001471','001333']
for dir in dirs:

     
    dataDir='/data1/IXPE/data/misureDU_2/'+dir
    print("ecco i files:")
    n_files=4
    if dir=='001461' or dir=='001471':
       n_files=8

    filenames=get_files(dataDir,n_files)
    print ("primi 4 files= ",filenames)
 
    out_dir=base_dir+dir+'/'


    n_iter=1

    zero_thr=bestZeroThr[dict_energy[dir]]
    moma1_thr=bestMoma12[dict_energy[dir]]
    moma2_thr=bestMoma12[dict_energy[dir]]
    dmin=bestDmin[dict_energy[dir]]
    weight_scale=bestws[dict_energy[dir]]

    print("dir=",dir," E = ",dict_energy[dir]," zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
    work_dir=out_dir+str(n_iter)
        
    #lanciaRecon_hct.submit_recon(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir)

    filename=work_dir+'/config_simo.txt'
    myLogfile=open(filename,'w')  
    myLogfile.write("zeroThreshold= "+str(zero_thr)+" moma1_thr= "+str(moma1_thr)+" moma2_thr= "+str(moma2_thr)+"  dmin= "+str(dmin)+" dmax= "+str(dmax)+" w_scale= "+str(weight_scale) )
    



    
    
    

