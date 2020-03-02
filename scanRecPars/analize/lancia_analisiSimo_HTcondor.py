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


# scan ZeroMoma12
#base_dir='/data1/maldera/IXPE_work/rec_optimization/scanZeroMoma12/'
#max_iter= 138+1


# scanRmin-Ws
base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/'
max_iter= 132+1





dirs=['001333','001361',  '001388',  '001416',  '001436',  '001461',  '001471']
#dirs=['001471']







stringa1='#!/bin/bash \n \nFILE="/var/log/All_done.txt"\nwhile true \ndo  \n   if [ -f $FILE ]; then \n      echo "File $FILE exists... let\'s go!" \n'
stringa2='\n      break  \n  \n   else  \n      echo "File $FILE does not exist...  waiting 30 secs" \n      sleep 30\n   fi  \ndone  \n'    




def run_command(cmd):
         
        print('going to run: ',cmd)
        subprocess.call(cmd,shell=True)



def crea_jobScript(command_string,work_dir):

         
         #command_string=command_string+' \n'
         nomeFile=work_dir+'/job_script_analisi.sh'
         miofile = open(nomeFile,'w')       # apre il file in scrittura
         miofile.write(stringa1)
         miofile.write(command_string) # !!!!!!!!!!!!!!!!!!!!!1 actual command!!!!!! 
         miofile.write(stringa2)
         miofile.close()               
         run_command('chmod a+x '+ nomeFile)


def lancia_jobHtc(work_dir, cmd_string):

        run_command('cp  /home/users/maldera/IXPE/ixpeTO/scanRecPars/analize/submit_analisi '+work_dir+'/.')
        crea_jobScript(cmd_string,work_dir)
        run_command('cd '+work_dir+' &&   condor_submit submit_analisi')
        

        





n_iter=1
#inizio scan su zero_thr
for dirName in dirs:

    out_dir=base_dir+dirName+'/'
    
    for  i in range (1,max_iter):
        work_dir=out_dir+str(i)
        files_string=work_dir+'/*.lv1_recon.fits'                
        cmd_string='      python /home/users/maldera/IXPE/ixpeTO/scanRecPars/analize/plotAll_simo.py '+files_string+' -o '+work_dir+'/'

        lancia_jobHtc(work_dir, cmd_string)        
                
        n_iter=n_iter+1

        