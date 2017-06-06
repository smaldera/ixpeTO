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
