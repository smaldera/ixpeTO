import subprocess

# per ogni ricostruzione, pernde tutti i files recon a tutte le enernge e lancia plotAll_PHA che divide in bin di L/W





# scanRmin-Ws
base_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/'
max_iter= 132+1



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
        

        




    
for  i in range (1,max_iter):
#for  i in range (1,2):

                

        #work_dir=out_dir+str(i)

        out_dir='/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/TLbins/rec_'+str(i)+'/'
        work_dir=out_dir
        run_command('mkdir -p '+work_dir)
        

        run_command('cp '+ base_dir+'*/'+str(i)+'/config_simo.txt '+out_dir+'.')
        
        

        
        files_string=base_dir+'*/'+str(i)+'/*_recon.fits'
        #/data1/maldera/IXPE_work/rec_optimization/scanRmin-Ws/*/30/*_recon.fits

        print ("file_string = ",files_string)
        
        cmd_string='      python /home/users/maldera/IXPE/ixpeTO/scanRecPars/analize/plotAll_simoLenght.py '+files_string+' -o '+out_dir

        print ("cmd=",cmd_string)
        
        lancia_jobHtc(work_dir, cmd_string)        
                        
