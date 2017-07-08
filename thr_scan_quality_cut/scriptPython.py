#python scriptPython.py /data1/IXPE/data/MC_namesOK/MDAT_files_namesOK/sim_7Kev.mdat -z 3



import os 



import argparse
formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument('infile', type=str,
                    help='the input binary file')
parser.add_argument('-z', '--zero-suppression', type=int, default=5,
                    help = 'zero-suppression threshold')


args = parser.parse_args()
FILE_PATH = args.infile
zero_sup = args.zero_suppression

print os.path.basename(FILE_PATH)
print 'ZERO_SUPPRESSION set to ', zero_sup, ' ADC'



folder_name = os.path.basename(FILE_PATH).replace('.mdat','')



for i in range (zero_sup, 21):
    for j in range (zero_sup, 21):
        cmd='/home/users/software/IXPE/gpdsw/bin/ixperecon --input-files ' + FILE_PATH + ' --output-folder /data1/IXPE/data/MC_namesOK/recon_FITS/' + folder_name +'/' + ' --output-suffix zero_%s' %zero_sup + '_th1_%s' %i + '_th2_%s' %j + ' --threshold %s' %zero_sup + ' --moma1-threshold %s' %i + ' --moma2-threshold %s' %j + ' -n 250000'  #+str(i)
        print "sto per eseguire ",cmd
        os.system(cmd)










'''
if __name__ == "__main__":


    import argparse
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter)
    parser.add_argument('infile', type=str,
                        help='the input binary file')
    parser.add_argument('-z', '--zero-suppression', type=int, default=5,
                        help = 'zero-suppression threshold')

    args = parser.parse_args()

    #FILE_PATH = os.path.join(os.environ['GPDSWROOT'], 'Recon', 'data', 'test_fe_500evts.mdat')
    FILE_PATH = args.infile

    f = TFile("thr_scan_%s_phi.root" %os.path.basename(FILE_PATH).replace('.mdat',''), "recreate") #NEW

    test(args.zero_suppression)

    f.Close()
'''


