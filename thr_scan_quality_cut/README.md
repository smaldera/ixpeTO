FOLDER: thr_scan_quality_cut

CONTENT USAGE

==================================
0. Environment settings
==================================
First of all, do: $source /path/to/ixpesw/setup.sh


==================================
1. MDAT files thr scan
==================================

--> $python scriptPython.py path/to/in_file.mdat -z zero_suppressioj_threshold

e.g.: $python scriptPython.py /data1/IXPE/data/MC_namesOK/MDAT_files_namesOK/sim_7Kev.mdat -z 3

Output: FITS files, one for each couple of thresholds
e.g.: sim_7keV_zero_3_th1_18_th2_9.fits


==================================
2. Eccentricity (M2TL) histograms
==================================

For each fits file (at a certain energy), plot TRK_M2T / TRK_M2L histograms

--> $python -W ignore M2TL_scan.py FITS_folder_name -z zero_suppression_threshold

e.g.: $python -W ignore M2TL_scan.py sim_7keV -z 3

N.B.: M2TL_scan.py must be in the folder in which the directories (one for each energy) containing FITS files are stored

Output: ROOT file, one for each energy, containing M2TL histogram for each couple of thresholds
e.g.: thr_scan_sim_6keV_M2TL.root


==================================
3. M2TL quality cut + TRK_PHI2 histograms
==================================

For each couple of thresholds and at a fixed energy, apply the quality cut on eccentricity (M2TL) and plot TRK_PHI2 histogram (using only the events which pass the quality cut)

--> $python set_quality_cut.py FITS_folder_name path/to/thr_scan_M2TR_file.root -z zero_suppression_threshold

e.g.: $python set_quality_cut.py sim_6keV thr_scan_sim_6keV_M2TL.root -z 3

N.B.: M2TL_scan.py must be in the folder in which the directories (one for each energy) containing FITS files are stored

Output: ROOT file, one for each energy, containing TRK_PHI2 histograms for each couple of thresholds
e.g.: thr_scan_sim_6keV_PHI2.root


==================================
4. Modulation factor analysis
==================================

Study the modulation factor for each of the TRK_PHI2 histogram in thr_scan_aaaa_PHI2.root

--> $python modulation_factor_analysis.py path/to/thr_scan_PHI2_file.root -z zero_suppression_threshold

e.g.: $python modulation_factor_analysis.py thr_scan_sim_PHI2_qcut/thr_scan_sim_9keV_PHI2.root -z 3

Output: ROOT file, one for each energy, containing (for each couple of thresholds):
    - phi distribution histograms FITTED
    - modulation factors scatter plot and 3D plot (1st and 2nd threshold on the x and y axes)
    - reduced chi square and probability histograms
    - modulation factor distribution histogram
e.g.: mod_fact_sim_6keV_PHI2.root






plot_modFactor.py

python script that reads the phi histograms created by the threshold scan and creates a plot for mu vs energy for the "standard" case (i.e. thresolds =5,5). It also plots the difference max_mu-std_mu vs energy.

the script needs 2 arguments: a config file and the name of the out root file. The config file is the list of the root files to be alayzed and its energy es:: histograms/mod_fact_root_files_SIM/mod_fact_sim_2_5keV_L_phi.root 2.5 ... histograms/mod_fact_root_files_SIM/mod_fact_sim_2keV_L_phi.root 2

USAGE: python -i plot_modFactor.py imputFileList.txt outFileName.root


