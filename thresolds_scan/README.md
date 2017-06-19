FOLDER: threshold_scan


CONTENT USAGE

--> first of all, do:
  $source /path/to/ixpesw/setup.sh
  
--> to read and analyse MDAT file:
  $python mom_analysis_threshold_scan_classe.py inFile.mdat
  
  The output file will be  -->  thr_scan_inFile_phi.root
  You'll find phi distribution histograms in it for each double threshold
  
--> to do the analysis of the modulation factors:
  $python modulation_factor_analysis_classe2.py thr_scan_inFile_phi.root
  
  The output file will be  -->  mod_fact_inFile_phi.root
  You'll find
      - phi distribution histograms FITTED
      - modulation factors scatter plot and 3D plot (1st and 2nd threshold on the x and y axes)
      - reduced chi square histogram

--> to study the projection on the two axes of the scatter plot, do:
  $python modulation_factor_analysis_matplotlib.py thr_scan_inFile_phi.root





--------------------------------------
plot_modFactor.py

python script that reads the phi histograms created by the threshold scan and creates
a plot for mu vs energy for the "standard" case (i.e. thresolds =5,5). It also plots
the difference max_mu-std_mu vs energy.

the script needs 2 arguments: a config file and the name of the out root file.
The config file is the list of the root files to be alayzed and its energy
es::
histograms/mod_fact_root_files_SIM/mod_fact_sim_2_5keV_L_phi.root   2.5
...
histograms/mod_fact_root_files_SIM/mod_fact_sim_2keV_L_phi.root     2



USAGE:  python -i  plot_modFactor.py  imputFileList.txt  outFileName.root

