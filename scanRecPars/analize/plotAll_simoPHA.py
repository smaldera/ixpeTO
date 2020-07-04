#***********************************************************************
# Copyright (C) 2017 the Imaging X-ray Polarimetry Explorer (IXPE) team.
#
# For the license terms see the file LICENSE, distributed along with this
# software.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#***********************************************************************
from __future__ import print_function, division


__description__ = 'analyze runs'

import matplotlib
matplotlib.use('Agg')  # questo fa funzionare matplotlib senza interfaccia grafica (es su un server... ) !!!!!!!! OKKIO!!!!

import numpy

import subprocess


from gpdswpy.binning import ixpeHistogram1d,ixpeHistogram2d
from gpdswpy.logging_ import logger
from gpdswpy.dqm import ixpeDqmTask, ixpeDqmArgumentParser
from gpdswpy.fitting import fit_gaussian_iterative, fit_histogram,  fit_modulation_curve
from gpdswpy.filtering import full_selection_cut_string, energy_cut, cut_logical_and
from gpdswpy.modeling import ixpeFe55
from gpdswpy.matplotlib_ import plt
from gpdswpy.tasks.pha_trend import pha_trend
from scipy.interpolate import InterpolatedUnivariateSpline


parser = ixpeDqmArgumentParser(description=__description__)


parser.add_argument('--normalize', action='store_true', default=False,
                    help='normalize the pulse height distribution area to 1')
parser.add_argument('--fit', action='store_true', default=False,
                    help='fit the pulse height distribution')
parser.add_argument('--fit-model', type=str, help='fit model',
                    default='gauss', choices=['gauss', 'Fe55'])
parser.add_argument('--fit-min', type=float, nargs='?', default=0.,
                    help='fit lower edge for Fe55 model')
parser.add_argument('--fit-max', type=float, nargs='?', default=40000.,
                    help='fit upper edge for Fe55 model')
parser.add_argument('--nsigma', type=float, nargs='?', default=1.5,
                    help='fitting range in number of sigma for gauss model')
parser.add_argument('--seconds-per-bin', type=float, nargs='?', default=600,
                    help='number of seconds in one temporal bin')
parser.add_argument('--correct-peak-drift', action='store_true',
                    default=False,
                    help='Correct spectrum for peak drift')

parser.add_argument('--degrees', action='store_true', default = False,
                    help='plot the modulation curve with angles in degrees')

parser.add_pha_options()
#parser.set_defaults(pha_expr='TRK_PI')
parser.set_defaults(pha_expr='TRK_PHA')

parser.add_cut_options()



### PARAMETRI....

cut_base='((abs(TRK_BARX) < 7.000) && (abs(TRK_BARY) < 7.000)) && ( (NUM_CLU > 0) && (LIVETIME > 15)  )'
ecut='(TRK_PI > 3700) && (TRK_PI<30000)'
track_size_cut='(TRK_SIZE > 0)'



pha_min =0
pha_max =50000
pha_bins =200 




def find_bins(run,  expr, cut):
    """
    """
    vals = run.values(expr, cut)
    deltabin = 0.0001
    # _nbinsqt = int(((max(vals) - min(vals)) / deltabin) + 1)
    _nbinsqt =int(((30 - min(vals)) / deltabin) + 1)
    

    nbinsqt = min([_nbinsqt, 5000000])
    logger.info('M2L/M2T histo n. bins %s' % nbinsqt)
    binning = numpy.linspace(min(vals), max(vals), nbinsqt)
    #binning = numpy.linspace(0, 2, 10)
   
     
    #hist = ixpeHistogram1d(binning, vals, xtitle=expr)
    hist = ixpeHistogram1d(binning, vals, xtitle=expr)

    cdf=hist.cdf()
    
    print("type of cdf=",type(cdf))

    my_bins=[0]
    N=200000
    k=1
    for x in binning:
       # print ("cdf(",x," )=",cdf(x)* hist.weights_sum)
        n=cdf(x)* hist.weights_sum
        if n>N*k:
            print ("k=",k,"x=",x," n=",n)
            my_bins.append(x)
            k+=1

    print ("my_bins=",my_bins)       




    hist2 = ixpeHistogram1d(numpy.array(my_bins), vals, xtitle=expr)
    hist2.plot(stat_box_position='upper right')
    plt.ylim(0,hist2.max_val()*1.1) #!!!!!!!!    


    for x in my_bins: 
       plt.axvline(x)

    return my_bins






class plotAll_simo(ixpeDqmTask):

    """
    """

    def do_run(self, **kwargs):
        """
        """


        #loop sui bin di L/W

        #binning = numpy.linspace(0, 2, 10)
        cut_base2=cut_logical_and(cut_base, ecut, track_size_cut )
        expr = 'TRK_M2L/TRK_M2T'
        #self.add_plot('moments ratioALL', hist, figure_name='moments ratioALL')
        my_bins=find_bins(self.run_list, expr, cut_base2)
        plt.savefig(kwargs.get('output_folder')+'/dist_ratioLW_ALL.png')
                                             
        cut_LW=''
        #for i in range(0,len(binning)-4):
        for i in range(0,len(my_bins)-1):
                     
            #cut_LW="( (  (TRK_M2L/TRK_M2T)>"+str(10**(binning[i]))+") && ( (TRK_M2L/TRK_M2T)<"+str(10**(binning[i+1]))+"))"
            cut_LW="( (  (TRK_M2L/TRK_M2T)>"+str(my_bins[i])+") && ( (TRK_M2L/TRK_M2T)<"+str(my_bins[i+1])+"))"

             
            print("cut_LW=",cut_LW)
        
            outFolder=kwargs.get('output_folder')+'/LWbin_'+str(i)+'/'

            #creo outFolder:
            cmd='mkdir -p '+outFolder
            print ("going to run: ",cmd)
            subprocess.call(cmd,shell=True)
              
        
            pha_expr ='TRK_PHA' #        kwargs.get('pha_expr')
            pha_binning = numpy.linspace(pha_min, pha_max, pha_bins)
            pha_title='Pulse height [ADC counts]'
        
            hist = ixpeHistogram1d(pha_binning, xtitle=pha_title)
           
                                        
            #cut=cut_base
            print ("base_cut=",cut_base)
            print ("cut_LW=",cut_LW)
         
            cut= cut_logical_and(cut_base2,cut_LW )
            print ("cut_final = ",cut)

            
        
            logger.info('Filling the histograms...')
            pha = self.run_list.values(pha_expr, cut) #!!!! qua fa il lavoro!!! 
            hist.fill(pha)                            
            print ("n_ev rimasti = ",len(pha))
            self.add_plot('pha_spectrum_'+str(i), hist, figure_name='pha_spectrum_'+str(i),  stat_box_position=None, label=kwargs.get('label'),  save=False)


   
            self.save_figure('pha_spectrum_'+str(i), overwrite=True)
            
            plt.savefig(outFolder+'PHA_spectrum1_'+str(i)+'.png')

            pha_mean=hist.mean[0]
            pha_rms=hist.rms[0]

            ##################################
            # plot dist L/W dopo i tagli:
            
            nBinsLW =int(((30 - 1) / 0.01) + 1)                        
            LW_binning=numpy.linspace(1, 30, nBinsLW)            
            histLW = ixpeHistogram1d(LW_binning, xtitle="L/W")
            LW = self.run_list.values("TRK_M2L/TRK_M2T", cut) #!!!! qua fa il lavoro!!!
            histLW.fill(LW)
            self.add_plot('LW_'+str(i), histLW, figure_name='LW_'+str(i), stat_box_position='upper right',  save=False)
            self.save_figure('LW_'+str(i), overwrite=True)
            plt.savefig(outFolder+'LW_'+str(i)+'.png')
            LW_mean=histLW.mean[0]
            LW_rms=histLW.rms[0]
            print("mean LW=",histLW.mean[0], "rms=",histLW.rms[0])                        
            
            


            ###################################3
            # mappa bary
                  
            x = self.run_list.values('TRK_BARX', cut)
            y = self.run_list.values('TRK_BARY', cut)
            x_min=-8
            x_max=8
            y_min=-8
            y_max=8
            nside=320
            x_edges = numpy.linspace(x_min, x_max, nside +1)
            y_edges = numpy.linspace(y_min, y_max, nside +1)
            hist_map = ixpeHistogram2d(x_edges, y_edges,  xtitle='x [mm]', ytitle='y [mm]')
            hist_map.fill(x, y)
            self.add_plot('bary map_'+str(i), hist_map, figure_name='bary_map_'+str(i))
            plt.savefig(outFolder+'bary_map_'+str(i)+'.png')
        
            ###################################3
            # istogramma ph1

            phi1= self.run_list.values('numpy.degrees(TRK_PHI1)', cut)
        
            edge=180
            #nbins=360
            nbins=180

            
            ang_binning = numpy.linspace(-edge, edge, nbins + 1)
            #print ("ang_bnins = ",ang_binning)
        
            hist_phi1 = ixpeHistogram1d(ang_binning, xtitle='deg')
            hist_phi1.fill(phi1)
            self.add_plot('hist_phi1_'+str(i), hist_phi1, figure_name='hist_phi1_'+str(i))
            # fitto phi1:
            fit_model1 = fit_modulation_curve(hist_phi1, xmin=-edge, xmax=edge, degrees=True,  verbose=kwargs.get('verbose'))
            self.add_plot('modulation_curve_'+str(i),  fit_model1,figure_name='hist_phi1_'+str(i), save=False, display_stat_box=kwargs.get('display_stat_box', True), position=kwargs.get('position', 'lower left'))



            phase1 = fit_model1.parameter_value('Phase')
            phase1_err = fit_model1.parameter_error('Phase')
            modulation1 = fit_model1.parameter_value('Modulation')
            modulation1_err = fit_model1.parameter_error('Modulation')
            chi2_1 = fit_model1.reduced_chisquare()
            plt.savefig(outFolder+'modulation_phi1_'+str(i)+'.png')
            ###################################3
            # istogramma ph2

            phi2= self.run_list.values('numpy.degrees(TRK_PHI2)', cut)
            #edge=180
            #nbins=360
        
            #ang_binning = numpy.linspace(-edge, edge, nbins + 1)
            hist_phi2 = ixpeHistogram1d(ang_binning, xtitle='deg')
            hist_phi2.fill(phi2)
            self.add_plot('hist_phi2_'+str(i), hist_phi2, figure_name='hist_phi2_'+str(i))
            #fit phi2
            fit_model2 = fit_modulation_curve(hist_phi2, xmin=-edge, xmax=edge, degrees=True,  verbose=kwargs.get('verbose'))
            self.add_plot('modulation_curve2_'+str(i),  fit_model2,figure_name='hist_phi2_'+str(i), save=False, display_stat_box=kwargs.get('display_stat_box', True), position=kwargs.get('position', 'lower left'))

            phase2 = fit_model2.parameter_value('Phase')
            phase2_err = fit_model2.parameter_error('Phase')
            modulation2 = fit_model2.parameter_value('Modulation')
            modulation2_err = fit_model2.parameter_error('Modulation')
            chi2_2 = fit_model2.reduced_chisquare()


            print("modulation2_err= ",modulation2_err)
            plt.savefig(outFolder+'modulation_phi2_'+str(i)+'.png')

            ################################################
            # count events:
            n_raw=self.run_list.num_events(cut_LW)
            n_final=self.run_list.num_events(cut)
            print ("n. raw events= ",n_raw)
            print ("n. final",n_final )
            
        


            #scrivi outfile

            
            nomefileout= outFolder+'prova_outLW.txt'
            print("nomefile_out= ",nomefileout )
            
            
            out_string=str(pha_mean)+' '+str(pha_rms)+' '+str(phase1)+' '+str(phase1_err)+' '+str(modulation1)+' '+str(modulation1_err)+' '+str(chi2_1)+' '+str(phase2)+' '+str(phase2_err)+' '+str(modulation2)+' '+str(modulation2_err)+' '+str(chi2_2)+' '+str(n_final)+' '+str(n_raw)+' '+str(my_bins[i])+' '+str(my_bins[i+1])+' '+str(LW_mean)+' '+str(LW_rms)
            
    
            with open(nomefileout, 'w') as miofile:
                miofile = open(nomefileout,'w')
                miofile.write(out_string)
                miofile.close() # !!!!! il file e' bufferizzato, e riempito solo alla chiusura (o chiamando file.flush). messo con with dovrebbe essere chiuso e scritto comunque quando esce dal loop
            
       
              
        
if __name__ == '__main__':
    args = parser.parse_args()
    opts = vars(args)
    opts['file_type'] = 'Lvl1a'
      
    task = plotAll_simo(*args.infiles, **opts)
    task.run(**opts)

    #if args.__dict__['show']:
    #plt.show()
