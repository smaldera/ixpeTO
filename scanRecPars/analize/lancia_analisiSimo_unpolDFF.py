# script per analizzare le ricistruzioni di ogni parametro
# versione per lo scan Unpol DFF (cambia solo l'aggiunta dei files dalle due direcoty eps)
#OKKIO  l'ouput vinene salvato nella prima dir eps  eps_dir[dirName][0]


import subprocess

#Rmin
base_dir='/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanDmin/'
max_iter=17+1
#max_iter=2                                                



#scanMoma1Thr
#base_dir='/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanMoma1/'
#max_iter=20+1

#scanMoma2Thr
#base_dir='/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanMoma2/'
#max_iter=20+1 

##exponential weight
#base_dir='/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scan_ws/'
#max_iter=21

# soglia clustering  
#base_dir='/data1/maldera/IXPE_work/rec_optimization/unpolDFF/scanZeroThr/'
#max_iter=19+1


#energy_dirs=['2.04','2.29', '2.70', '2.98', '3.69', '5.89']
energy_dirs=['5.89']
eps_dir={'2.04':['000647','000658']  ,'2.29':['000669','000677'], '2.70':['000686','000691'], '2.98':['000704', '000726'], '3.69':['000733','000740'], '5.89':['000744','000752']}   # dir angolo1 e angolo2






n_iter=1
#inizio scan su zero_thr
for dirName in energy_dirs:

    out_dir=base_dir+dirName+'/'

    
    for  i in range (1,max_iter):
    
      
        work_dir1=out_dir+'/'+eps_dir[dirName][0]+'/'+str(i)
        work_dir2=out_dir+'/'+eps_dir[dirName][1]+'/'+str(i)
        
        print("dir_1=",work_dir1)
        print("dir_2=",work_dir2)

        
        
        files_string=work_dir1+'/*.lv1_recon.fits '+work_dir2+'/*.lv1_recon.fits'                
        cmd='python plotAll_simo.py '+files_string+' -o '+work_dir1+'/'
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)

        
        
        n_iter=n_iter+1
