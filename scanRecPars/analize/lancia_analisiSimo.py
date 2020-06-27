import subprocess


#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanZeroThr_v2/'
#max_iter=19+1



#scanMoma1Thr
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanMoma1Thr/'
#max_iter=21

#scanMoma2Thr
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanMoma2Thr/'
#max_iter=21                                                                                                                                      

#Rmin
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin/'
#max_iter=18

#Rmax
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmax/'
#max_iter=22


#exponential weight
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanWs/'
#max_iter=21



# scan 3D  Dmin-ws + moma1,2
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/001333/'
base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/001361/'
max_iter=48+1





#dirs=['001333','001361',  '001388',  '001416',  '001436',  '001461',  '001471']
#dirs=['001333']


dirs=['moma1_20','moma1_22','moma1_24','moma1_26','moma1_28','moma1_30']








n_iter=1
#inizio scan su zero_thr
for dirName in dirs:

    out_dir=base_dir+dirName+'/'
    #out_dir=base_dir # no loop!!!!!!!

    
    
    for  i in range (1,max_iter):
    
      
        work_dir=out_dir+str(i)
        files_string=work_dir+'/*.lv1_recon.fits'                
        cmd='python plotAll_simo.py '+files_string+' -o '+work_dir+'/'
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)

        
        
        n_iter=n_iter+1
