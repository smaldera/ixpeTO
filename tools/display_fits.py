#!/bin/sh
import os
import argparse
from astropy.table import Table
import numpy as np
import pandas
import matplotlib.pyplot as plt
import logging

""" Configure the main terminal logger.
    [To be moved to another config file to import...]
"""
logger = logging.getLogger('ixpetool')
logger.setLevel(logging.DEBUG)
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.DEBUG)
consoleFormatter = logging.Formatter(">>> %(message)s")
consoleHandler.setFormatter(consoleFormatter)
logger.addHandler(consoleHandler)



__description__ = 'Producing some useful plots from a fits file'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-in', '--infile', type=str, required=True,
                    help='the input .fits file')
PARSER.add_argument('-v', '--vars', type=str, nargs='+', required=False,
                    default=None,
                    help='names of variables of which you want to display the histogram'+\
                        'example: python display_fits.py -v "TRK_SKEW" "TRK_PHI2" "TIME"')
PARSER.add_argument('-s', '--scatter', type=str, nargs='+', required=False,
                    default=None,
                    help='names of the 2 variables of which you want to display the scatter plot'+\
                        'example: python display_fits.py -s "TRK_SKEW" "TRK_PHI2"')


def display_variable_histo(pandas_dataframe, variable_name):
    """Simple function to display a variable histogram
    """
    plt.figure(facecolor='white')
    pandas_dataframe[variable_name].hist(bins=100)
    plt.ylabel('Counts')
    plt.xlabel(variable_name)

def display_scatter_plot(pandas_dataframe, variable_name1, variable_name2):
    """Simple function to display the scatter plot of 2 variables
    """
    plt.figure(facecolor='white')
    plt.hist2d(pandas_dataframe[variable_name1], pandas_dataframe[variable_name2],
                bins=50)
    plt.xlabel(variable_name1)
    plt.ylabel(variable_name2)
    cb = plt.colorbar()
    cb.set_label('counts')
    plt.grid()

def main(**kwargs):
    """Main function:
            1) get FITS file and turn it into a pandas DataFrame:
            2) if: -v var1, ..., varN, then histograms of var1, ... 
               up to varN are shown.
            3) if: -s var1, var2, then the scatter plot between var1
               and var2 is shown.
    """
    assert(kwargs['infile'].endswith('.fits'))
    data_fits = Table.read(kwargs['infile'], format='fits')
    df = data_fits.to_pandas()
    logger.info('Now file %s is a pandas DataFrame!' %os.path.basename(kwargs['infile']))
    logger.info('List of Columns:\n %s\n'%str(np.array(list(df.columns))))
    if kwargs['vars'] is not None:
        for col in kwargs['vars']:
            logger.info('Display variable %s'%col)
            display_variable_histo(df, col)
    if kwargs['scatter'] is not None:
        var1 = kwargs['scatter'][0]
        var2 = kwargs['scatter'][1]
        logger.info('Plotting %s vs %s'%(var1, var2))
        display_scatter_plot(df, var1, var2)
    if kwargs['scatter'] is not None or kwargs['vars'] is not None:
        plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    main(**args.__dict__)
