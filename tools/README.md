FOLDER: tools


PROGRAMS DESCRIPTION:

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

USAGE: $python -i read_fits.py path/to/FITSfile.fits output_file_name.root


===========================================================
NAME: read_TTree_hist_ASKED.py

--> script per leggere il TTree contenuto nel file ROOT (che
    deve essere prodotto usando il programma read_fits.py)
    e per disegnare (con ROOT) l'istogramma richiesto

USAGE: $python -i read_TTree_hist_ASKED.py file_name.root requested_plot



