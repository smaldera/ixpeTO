
import os
import numpy as np
import time
import subprocess
import time

def eff(E):

   return 2.7**(-0.63*E--0.071) 


N0=1000.
#N0=100.

n_iter=4
energy=[4.5]
#energy=[3]

angles=[-180]


#energy=[3]
#angles=[-175]



for nSims in range (0,1):

   for ene in energy :

      n=int(N0/eff(ene))
      print "n eventi da simualre =",n

      for angle in angles:
         print "energy = ",ene," angolo = ",angle
         n_iter=n_iter+1
        
         #seed=int(round(time.time() ) )+n_iter
         seed=int(repr(time.time()).split('.')[1])+n_iter

         
         miofile = open('job_script.sh','w')       # apre il file in scrittura
         miofile.write('#!/bin/bash  \n\n')

         miofile.write('FILE="/var/log/All_done.txt"   \n')
    
         miofile.write('while true  \n')
         miofile.write('do  \n')
         miofile.write('  if [ -f $FILE ]; then   \n')
         miofile.write('        echo "File $FILE exists."  \n')
    
         #miofile.write( " /home/users/maldera/IXPE/gpdsw/bin/ixpesim  --num-events "+str(n)+" --src-spectrum line --src-energy "+str(ene)+" --src-polarized 1  --src-pol-angle  90   --src-theta "+ str(angle)+"  --src-pos-x 6 --output-file  sim_"+str(ene)+"KeV_"+str(angle)+"deg.fits --output-root-file sim_"+str(ene)+"KeV_"+str(angle)+"deg.root  --log-file sim_"+str(ene)+"KeV_"+str(angle)+"deg.root   \n") # write scrive solo stringhe!!! 
         
         miofile.write( "/home/users/maldera/IXPE/gpdsw/bin/ixpesim  --num-events "+str(n)+" --src-spectrum line --src-energy "+str(ene)+" --src-polarized 1    --src-theta "+ str(angle)+"  --src-pos-x 0 --output-file  sim_"+str(ene)+"KeV_"+str(angle)+"deg.fits  --log-file sim_"+str(ene)+"KeV_"+str(angle)+"deg.log  --random-seed "+str(seed) + " --src-pol-angle 45   --src-sigma 4  --src-pos-z 200    \n") # write scrive solo stringhe!!! 
        

         miofile.write('        break  \n')
         miofile.write('  else  \n')
         miofile.write('        echo "File $FILE does not exist... waiting 30 secs" \n')
         miofile.write('        sleep 30   \n')
         miofile.write('  fi  \n')
         miofile.write('done  \n')
 

         miofile.close()


         myBash=open('myBash.sh','w')
         myBash.write('#!/bin/bash  \n\n')
         myBash.write('mkdir -p  '+str(n_iter)+" \n") 
         myBash.write('cp  submitJob_simo '+str(n_iter)+"/.  \n") 
         myBash.write("chmod a+x  job_script.sh \n" )
         myBash.write("mv  job_script.sh  "+str(n_iter)+"/. \n" )
    
         myBash.write(" cd  "+str(n_iter)+" \n" )
         myBash.write(" echo 'going to run... pwd' \n" )
         myBash.write("pwd \n" )
         
         myBash.write("condor_submit submitJob_simo \n" )  
         myBash.write(" cd  .. \n" )
         myBash.close()
   
       

         subprocess.Popen(['/bin/bash', '-c', 'source myBash.sh'])
         print "fatto...  "
         time.sleep(1)
 


    
print "the end... "

