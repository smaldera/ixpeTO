
import os
import numpy as np
import time
import subprocess
import time

stringa1='#!/bin/bash \n \nFILE="/var/log/All_done.txt"\nwhile true \ndo  \n   if [ -f $FILE ]; then \n      echo "File $FILE exists... let\'s go!" \n'
stringa2='      break  \n  \n   else  \n      echo "File $FILE does not exist...  waiting 30 secs" \n      sleep 30\n   fi  \ndone  \n'    






def crea_cmd(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev):


   mystring='ixperecon '+filenames+'   --output-folder '+output_folder+' --threshold '+str(zero_thr)+' --moma1-threshold '+str(moma1_thr)+' --moma2-threshold '+str(moma2_thr)+' --coherent-noise-offset '+str(coer_noise_offset)+'  --trigger-minicluster-offset '+str(trig_minicluster_offset)+'  --s '+suffix+'  --min-track-hits '+str(min_track_hits)+'   --min-density-points '+str(min_densy_points)+'  --dmin '+str(dmin)+'  --dmax '+str(dmax)+'  --weight-scale '+str(weight_scale)+' --truncate-lsb '+str(truncate_lsb)+'   --log-file '+logfile+' -n '+str(n_max)+ ' -f '+str(first_ev)+' \n'    


   print ('cmd= ',mystring)
   return mystring



def crea_jobScript(command_string):

         
         #command_string=command_string+' \n'
         miofile = open('job_script.sh','w')       # apre il file in scrittura
         miofile.write(stringa1)
         miofile.write(command_string) # !!!!!!!!!!!!!!!!!!!!!1 actual command!!!!!! 
         miofile.write(stringa2)
         miofile.close()               



def lancia_condorJob(outDir):

         cmd='mkdir -p  '+outDir+' \n'
         print("going to run: ",cmd[:-1])
         subprocess.call(cmd,shell=True)

         cmd='cp  submitJob_simo '+outDir+'/.  \n'
         print("going to run: ",cmd[:-1])
         subprocess.call(cmd,shell=True)

         cmd='chmod a+x  job_script.sh   \n'
         print("going to run: ",cmd[:-1])
         subprocess.call(cmd,shell=True)

         cmd='mv  job_script.sh  '+outDir+'/. \n'
         print("going to run: ",cmd[:-1])
         subprocess.call(cmd,shell=True)
         
         cmd='cd '+outDir+' && condor_submit submitJob_simo  \n'
         print("going to run: ",cmd[:-1])
         subprocess.call(cmd,shell=True)
     
 


def submit_recon(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev,work_dir):

   """
   filenames="/data1/IXPE/data/misureDU_2/001333/ixpe2019100*_1118.lv1.fits"
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
   n_max=100
   first_ev=0
   work_dir='/data1/maldera/IXPE_work/rec_optimization/test1'
   """

   cmd_str= crea_cmd(filenames,output_folder,zero_thr,moma1_thr,moma2_thr,coer_noise_offset,trig_minicluster_offset,suffix,min_track_hits,min_densy_points,dmin,dmax,weight_scale, truncate_lsb, logfile,n_max,first_ev)
   crea_jobScript(cmd_str)
   lancia_condorJob(work_dir)


   print ("job sottomesso!")



   
def submit_sim(num_events, src_spectrum, src_energy, src_polarized,  src_theta, src_pos_x, src_pos_y, output_file, log_file,  random_seed,  src_pol_angle,  src_sigma,   src_pos_z, work_dir): 

      
   cmd_ixpesim= 'ixpesim  --num-events '+str(num_events)+' --src-spectrum '+src_spectrum+' --src-energy '+str(src_energy)+" --src-polarized "+str(src_polarized)+ ' --src-theta '+ str(src_theta)+'  --src-pos-x '+str(src_pos_x)+' --output-file '+ output_file+ ' --log-file  '+ log_file +' --random-seed '+str(random_seed) + ' --src-pol-angle '+str(src_pol_angle)+' --src-sigma '+str(src_sigma)+' --src-pos-z '+str(src_pos_z)+' \n' 


   
   crea_jobScript(cmd_ixpesim)
   lancia_condorJob(work_dir)


   print ("job ixpesim  sottomesso!")

