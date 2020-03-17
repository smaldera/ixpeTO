import time
import lanciaRecon_hct


def eff(E):

   return 2.7**(-0.63*E--0.071) 


N0=400
#N0=100.

n_iter=0
#energy=[2,3,4,5,6,7,8]
energy=[7,8]       

angles=[-180]


#energy=[3]
#angles=[-175]

base_dir='/data1/maldera/IXPE_work/rec_optimization/simData/'


src_spectrum='line'
src_polarized=1
src_pos_x=0
src_pos_y=0
src_pol_angle=75
src_sigma=5
src_pos_z=400
#output_file='out_sim.fits'
#log_file='out_sim.log'


for nSims in range (1,30):

   for ene in energy :

      n=int(N0/eff(ene))
      print ("n eventi da simualre =",n)

      for angle in angles:
         print ("energy = ",ene," angolo = ",angle)
         n_iter=n_iter+1


         num_events=n
         src_energy=ene
         src_theta=angle
         
         #seed=int(round(time.time() ) )+n_iter
         random_seed=int(repr(time.time()).split('.')[1])+n_iter

         work_dir=base_dir+str(ene)+'KeV/angle'+str(angle)+'deg/'+str(nSims)+'/'
         output_file=work_dir+'out_sim.fits'
         log_file=work_dir+'out_sim.log'
         
         
#         miofile.write( "/home/users/maldera/IXPE/gpdsw/bin/ixpesim  --num-events "+str(n)+" --src-spectrum line --src-energy "+str(ene)+" --src-polarized 1    --src-theta "+ str(angle)+"  --src-pos-x 0 --output-file  sim_"+str(ene)+"KeV_"+str(angle)+"deg.fits  --log-file sim_"+str(ene)+"KeV_"+str(angle)+"deg.log  --random-seed "+str(seed) + " --src-pol-angle 45   --src-sigma 5  --src-pos-z 400    \n") # write scrive solo stringhe!!! 
        
 
         lanciaRecon_hct.submit_sim(num_events, src_spectrum, src_energy, src_polarized,  src_theta, src_pos_x, src_pos_y, output_file, log_file,  random_seed,  src_pol_angle,  src_sigma,   src_pos_z, work_dir)



         filename=work_dir+'/config_ixpeSim_simo.txt'
         myLogfile=open(filename,'w')
         myLogfile.write(str(num_events)+'  '+str(src_spectrum)+'  '+str(src_energy)+'  '+str(src_polarized)+'  '+str(src_theta)+'  '+str(src_pos_x)+'  '+str(src_pos_y)+'  '+str(output_file)+'  '+str(log_file)+'  '+str(random_seed)+'  '+str(src_pol_angle)+'  '+str(src_sigma)+'  '+str(src_pos_z)+ '\n' )
         myLogfile.close()
         
         
print ("the end... ")

