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


__description__ = 'PHA spectrum'

import matplotlib
#matplotlib.use('Agg')  # questo fa funzionare matplotlib senza interfaccia grafica (es su un server... )

import numpy

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
cut_sigma=3
cut_base='((abs(TRK_BARX) < 7.000) && (abs(TRK_BARY) < 7.000)) && ( (NUM_CLU > 0) && (LIVETIME > 15)  )'
ecut='(TRK_PI > 3700) && (TRK_PI<30000)'
track_size_cut='(TRK_SIZE > 0)'



pha_min =0
pha_max =50000
pha_bins =200 


def find_quantile(run, quantile, expr, cut):
    """
    """
    vals = run.values(expr, cut)
    deltabin = 0.01
    # _nbinsqt = int(((max(vals) - min(vals)) / deltabin) + 1)
    _nbinsqt =int(((30 - min(vals)) / deltabin) + 1)
    

    nbinsqt = min([_nbinsqt, 50000])
    logger.info('M2L/M2T histo n. bins %s' % nbinsqt)
    #binning = numpy.linspace(min(vals), max(vals), nbinsqt)
    binning = numpy.linspace(0, 2, 10)
   
     
    #hist = ixpeHistogram1d(binning, vals, xtitle=expr)
    hist = ixpeHistogram1d(binning, numpy.log10(vals), xtitle=expr)

     
#    q = hist.quantile(1. - quantile)
    hist.plot(stat_box_position='upper right')
    plt.yscale('log')

    
    #plot_xmax = hist.quantile(0.99)
    #plt.xlim(binning[0], plot_xmax)
    plt.ylim(0,hist.max_val()*1.1) #!!!!!!!!    
    #plt.ylim(0,30) #!!!!!!!!
    #plt.yscale('log')
    
    # plt.axvline(q)
   # return q



def peak_cut(model):
    """
    """
    peak = model.parameter_value('Peak')
    sigma = model.parameter_value('Sigma')
    #num_cut_sigma = kwargs.get('cut_sigma')
    num_cut_sigma = cut_sigma
    
    plt.axvline(peak - num_cut_sigma * sigma, label='selection region')
    plt.axvline(peak + num_cut_sigma * sigma)
    return energy_cut(peak, sigma, nsigma=num_cut_sigma)






class plotAll_simo(ixpeDqmTask):

    """
    """

    def do_run(self, **kwargs):
        """
        """


        #qua ci metto il loop sui bin di L/W

        binning = numpy.linspace(0, 2, 10)

        cut_LW=''
        for i in range(0,len(binning)-4):
            cut_LW="( (  (TRK_M2L/TRK_M2T)>"+str(10**(binning[i]))+") && ( (TRK_M2L/TRK_M2T)<"+str(10**(binning[i+1]))+"))"
            print("cut_LW=",cut_LW)
        
        
        
        
            pha_expr ='TRK_PHA' #        kwargs.get('pha_expr')
            pha_binning = numpy.linspace(pha_min, pha_max, pha_bins)
            pha_title='Pulse height [ADC counts]'
        
            hist = ixpeHistogram1d(pha_binning, xtitle=pha_title)
           
                                        
            #cut=cut_base
            print ("base_cut=",cut_base)
            print ("cut_LW=",cut_LW)
         
            cut_base2= cut_logical_and(cut_base,cut_LW, ecut, track_size_cut )
            logger.info('Full selection cut (for pha_spectrum) : %s' % cut_base2)
       
        
            logger.info('Filling the histograms...')
            pha = self.run_list.values(pha_expr, cut_base2) #!!!! qua fa il lavoro!!! 
            hist.fill(pha)                            
            print ("n_ev rimasti = ",len(pha))
            self.add_plot('pha_spectrum_'+str(i), hist, figure_name='pha_spectrum_'+str(i),  stat_box_position=None, label=kwargs.get('label'),  save=False)


   
            self.save_figure('pha_spectrum_'+str(i), overwrite=True)
            
            plt.savefig(kwargs.get('output_folder')+'PHA_spectrum1_'+str(i)+'.png')
        
            ###################################3
            # mappa bary e mappa punto impatto
            # mi serve un histo 2D... 

            # track_size_cut='(TRK_SIZE > 0)'
            cut2= cut_base2


            n_physical=self.run_list.num_events(cut_base2)
            n_ecut=self.run_list.num_events(cut2)
        
            ecut_efficiency = n_ecut/n_physical
            quantile = min(0.8/ecut_efficiency, 1.)
            expr = 'TRK_M2L/TRK_M2T'
        
            self.add_plot('moments ratio', hist, figure_name='moments ratio')
            min_mom_ratio = find_quantile(self.run_list, quantile, expr, cut2)
            plt.savefig(kwargs.get('output_folder')+'dist_ratioLW_'+str(i)+'.png')
           
            cut_final=cut2

        
            print ("cut_final = ",cut_final)
       
            x = self.run_list.values('TRK_BARX', cut_final)
            y = self.run_list.values('TRK_BARY', cut_final)
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
            plt.savefig(kwargs.get('output_folder')+'bary_map_'+str(i)+'.png')
        
            ###################################3
            # istogramma ph1

            phi1= self.run_list.values('numpy.degrees(TRK_PHI1)', cut_final)
        
            edge=180
            nbins=360
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
            plt.savefig(kwargs.get('output_folder')+'modulation_phi1_'+str(i)+'.png')
            ###################################3
            # istogramma ph2

            phi2= self.run_list.values('numpy.degrees(TRK_PHI2)', cut_final)
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
            plt.savefig(kwargs.get('output_folder')+'modulation_phi2_'+str(i)+'.png')
        
            ##############################################
            # rifaccio istogramma pha con tagli finali per avere la risuluzione!!!

            pha_binning = numpy.linspace(pha_min, pha_max, pha_bins + 1)
            
            hist = ixpeHistogram1d(pha_binning, xtitle=pha_title)
            #if (pha_expr == 'TRK_PI'):
            pha_title = 'PHA [ ADC counts]'
            hist2 = ixpeHistogram1d(pha_binning, xtitle=pha_title)
            pha2 = self.run_list.values(pha_expr, cut_final) #!!!! qua fa il lavoro!!! 
            hist2.fill(pha2)                            
            self.add_plot('pha_spectrum2_'+str(i), hist2, figure_name='pha_spectrum2',  stat_box_position=None, label=kwargs.get('label'),  save=False)

            #if kwargs.get('fit'):
            # if kwargs.get('fit_model') == 'gauss':
            nsigma = kwargs.get('nsigma')
        
            gauss_model2 = fit_gaussian_iterative(hist2, verbose=kwargs.get('verbose'), xmin=kwargs.get('fit_min'),  xmax=kwargs.get('fit_max'), num_sigma_left=nsigma,  num_sigma_right=nsigma) # n. iterazioni??
            
            self.add_plot('pha_spectrum_fit2_'+str(2),gauss_model2  , figure_name='pha_spectrum2',    save=True,       display_stat_box=kwargs.get('display_stat_box', True),    position=kwargs.get('position', 'upper left'))
       

            plt.savefig(kwargs.get('output_folder')+'PHA_spectrum2_'+str(i)+'.png')  # non riesco a salvare i plot usando i metodi della classe...  cosi' va...


        
            peak2 = gauss_model2.parameter_value('Peak')
            peak2_err = gauss_model2.parameter_error('Peak')
            resolution2 = gauss_model2.resolution()
            resolution2_err = gauss_model2.resolution_error()
            print("peak2 = ",peak2," +- ",peak2_err," res2 fwhm =",resolution2," +- ",resolution2_err)


        

        


            ################################################
            # count events:
            n_raw=self.run_list.num_events()
            n_physical=self.run_list.num_events(cut_base2)
            n_ecut=self.run_list.num_events(cut2)
            n_final=self.run_list.num_events(cut_final)
        
            print ("n. raw events= ",n_raw)
            print ("n. physical events (bary + livetime+ NUM_CLU ) = ", n_physical)
            print ("n. ecut (physical+pha_spectrum+_tkr_size) ",n_ecut )
            print ("n. final (physical+pha_spectrum+_tkr_size+axis ratio) ",n_final )
        
            print("eff_ecut= n_ecut/physical",float(n_ecut)/float(n_physical))
        
            print("eff= final/physical",float(n_final)/float(n_physical))
        


            #scrivi outfile

            print("out dir", kwargs.get('output_folder'))
            nomefileout= kwargs.get('output_folder')+'prova_out.txt'
            print("nomefile_out ",nomefileout )
            
            """
            out_string=str(peak2)+' '+str(peak2_err)+' '+str(resolution2)+' '+str(resolution2_err)+' '+str(phase1)+' '+str(phase1_err)+' '+str(modulation1)+' '+str(modulation1_err)+' '+str(chi2_1)+' '+str(phase2)+' '+str(phase2_err)+' '+str(modulation2)+' '+str(modulation2_err)+' '+str(chi2_2)+' '+str(n_raw)+' '+str(n_physical)+'  '+str(n_ecut)+'  '+str(n_final) 
    
            with open(nomefileout, 'w') as miofile:
                miofile = open(nomefileout,'w')
                miofile.write(out_string)
                miofile.close() # !!!!! il file e' bufferizzato, e riempito solo alla chiusura (o chiamando file.flush). messo con with dovrebbe essere chiuso e scritto comunque quando esce dal loop
            """
       
              
        
if __name__ == '__main__':
    args = parser.parse_args()
    opts = vars(args)
    opts['file_type'] = 'Lvl1a'
      
    task = plotAll_simo(*args.infiles, **opts)
    task.run(**opts)

    #if args.__dict__['show']:
    #plt.show()
