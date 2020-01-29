import subprocess



base_dir='/data1/maldera/IXPE_work/rec_optimization/scanZeroThr/'
dirs=['001361',  '001388',  '001416',  '001436',  '001461',  '001471']
#dirs=['001333']



n_iter=1
#inizio scan su zero_thr
for dir in dirs:

    out_dir=base_dir+dir+'/'

    
    for  i in range (1,19):
        work_dir=out_dir+str(i)
        files_string=work_dir+'/*.lv1_recon.fits' 
        cmd='python plotAll_simo.py '+files_string+' -o '+work_dir+'/'
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)

        
        
        n_iter=n_iter+1