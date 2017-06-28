suff relared to Convoluted Neural Network reonstruction


- sim2ascii.py:

  Reads the output of ixpesim (using the pyton wrapper of gpdsw ) and write the track on a txt file.
  For now it uses the mdat file to read the pixel info and the root file for the MC part (when the wrapper for the fits file will be avilable  we can swich to fits).
  The scrips call the clustering algorithm of ixperecon and writes only the pixle of the main track. 
  For each envent it writes:

  event_id n_pixels xConv yConv photonEnergy(MeV) phi_photoelectron  theta_photoelectron
  x_1 y_1  adc_1
  ...
  x_n y_n adc_n

  usage:
  source gpdsw/setup.sh  (you need to compile swig)
  python   sim2ascii.py  (input and output filenames are hardcoded)
  