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



