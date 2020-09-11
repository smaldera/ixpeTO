# lancia plotAll_simo.py sulle ricostruzioni bestparameters
# (generate con runBestParameters.py)



import subprocess



base_dir='/data1/maldera/IXPE_work/rec_optimization/bestParams/'
#iters=['best3d','std','para','ixperecon_para']
#iters=['ixpereconPara_LW','std']
#iters=['ixpereconPara_LxE']
#iters=['ixpereconPara_LW_new1']                                                                                                                 
#iters=['ixpereconPara_LW_new2']


#iters=['ixprecon_new']
#iters=['my_ixperecon_para']
#iters=['oldStd_ixperecon_para']
#iters=['my_ixperecon_para_v2'] 
#iters=['my_ixperecon_para_v3']
#iters=['my_ixperecon_para_v4']
iters=['my_ixperecon_para_v5']


dirs=['001333','001361',  '001388',  '001416',  '001436',  '001461',  '001471']





n_iter=1
#inizio scan su zero_thr
for dirName in dirs:

    out_dir=base_dir+dirName+'/'
        
    
    for  params in iters:
    
      
        work_dir=out_dir+params
        files_string=work_dir+'/*.lv1_reco*.fits'                
        cmd='python plotAll_simo.py '+files_string+' -o '+work_dir+'/'
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)

        
        
        n_iter=n_iter+1
