FOLDER: tools


PROGRAMS DESCRIPTION:

===========================================================
NAME: display_fits.py

--> script per leggere i file .fits (ixperecon output files)
e disegnare (con matplotlib) gli istogrammi delle variabili 
e scatter plot tra due variabili

USAGE: $python progName.py -in fileName.fits -v 'var1' 'var2'
        -s 'var1' 'var2'

EXAMPLE: 
    python tools/display_fits.py -in tools/data/xpol_2735_recon.fits -s 'TRK_SKEW' 'TRK_PHI2' -v 'TRK_SKEW' 'TRK_PHI2' 'TIME'


===========================================================
NAME: read_fits_hist.py

--> script per leggere i file .fits (ixperecon output files)
    e disegnare (con matplotlib) alcuni istogrammi

USAGE: $python progName.py fileName.fits requested_plot


===========================================================
NAME: read_fits.py

--> script per leggere i file .fits (ixperecon output files)
    e creare un ROOT.TTree i cui branch rappresentano le
    grandezze che entrano in gioco nella ricostruzione della
    traccia e nei quali sono salvati i dati relativi a tutti
    gli eventi contenuti nel file FITS dato in input.

USAGE: $python read_fits.py path/to/FITSfile.fits output_file_name.root


===========================================================
NAME: read_TTree_hist_ASKED.py

--> script per leggere il TTree contenuto nel file ROOT (che
    deve essere prodotto usando il programma read_fits.py)
    e per disegnare (con ROOT) l'istogramma richiesto

USAGE: $python read_TTree_hist_ASKED.py file_name.root requested_plot


===========================================================
NAME: read_TTree_hist.py

--> script per leggere il TTree contenuto nel file ROOT (che
    deve essere prodotto usando il programma read_fits.py)
    e per disegnare (con ROOT) quattro istogrammi hardcoded
    (non selezionabili da linea di comando)
--> at the end of the program, you can find (commented!!) an
    alternative much simpler basic code to plot four default
    hardcoded distribution histograms

USAGE: $python read_TTree_hist.py file_name.root



===========================================================
NAME: ixpeRecon_fits2rootSimo2.py

--> traduce un file fits creato da ixpeRecon in un TTree di root.
    il tree vine sui chiama "reconT" e viene salvato nel file di output specificato
   
USAGE: $python   ixpeRecon_fits2rootSimo2.py  file_recon.fits  outRootFile.root



===========================================================
NAME :mod_factors_fromRootSimo.py

--> small script to plot the modulation factor vs energy.
   imput files, energies and name of output files are hardcoded!!!  


USAGE: $python  mod_factors_fromRootSimo.py





===========================================================
NAME: readSimTree_Simo.cc

--> reads the root ttree resulting from the ixpe simulation.
    For each events draws 3 views of the energy deposit and a xy
    view of the different tracks. A canvas is saved for each event
    and it is saved in the out root file.
    
--> compile with:
    c++ -Wl,--no-as-needed   readSimTree_Simo.cc -o readSimTree_Simo.exe `root-config --libs --cflags`
   
USAGE:  ./readSimTree_Simo.exe inputfile.root  outfile.root  n_max_events (optional, default=100, <0 all events )


===========================================================
NAME: drawIxpeReconFits_simo.py

--> reads and plots recon fits file.
    Si possono aggiungere dei tagli (vedi sezione CUTS nel codice)


USAGE:  python drawIxpeReconFits_simo.py  recinFits_filesList.txt  outRootFile.root


