# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains some global metadata used throughout the brutus code to fit IFU data 
 cubes. Includes reference wavelengths, the speed of light, etc ...
 
 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import os

import matplotlib.pyplot as plt

# ---| Some basic parameters |------------------------------------------------------------

__version__ = '0.3.1'

# Where are we located ?
# Get the project directory from this file location !
brutus_dir = os.path.dirname(__file__)

# And get all the different useful locations associated
refdata_dir = os.path.join(brutus_dir,'..','reference_data')

# ---| Some constants |-------------------------------------------------------------------

c = 299792.458 # speed of light in km/s

# Some emission line rest-wavelengths (in air and Ansgtroem)
o2l = 3726.09
o2h = 3728.83
hd = 4101.73
hc = 4340.47
hb = 4861.32
o3l = 4958.91
o3h = 5006.73
he1 = 5875.66
nal = 5889.95
nah = 5895.92
o1 = 6300.304
n2l = 6548.04
ha = 6562.81
n2h = 6583.46
s2l = 6716.44
s2h = 6730.81

# And some names for the emission lines, useful for dealing with pyqz
elines_pyqz = {'Hb':[hb],
              '[OIII]':[o3h],
              'Ha':[ha],
              '[NII]':[n2h],
              '[SII]+':[s2h,s2l],
             }


# ---| Information about he file format for the different instruments |-------------------

ffmt = {'MUSE': {'data':1, 'error':2, 'badpix':None},
        #'WiFeS':{'data':0, 'error':1, 'badpix':2},
       }

# ---| Information about the different stellar models for ppxf|---------------------------

sl_models ={'MILES_ppxf_default':{'sl_loc':os.path.join(refdata_dir,
                                                        'stellar_templates',
                                                        'miles_models'),
                                  'nAges':26, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(1),np.log10(17.7828),26),
                                  'metal_str':['m1.71','m1.31','m0.71','m0.40','p0.00',
                                               'p0.22'],
                                  'fwhm':2.51,
                                }, 
            'MILES_Padova00_kb_1.30':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30'),
                                  'nAges':50, 'nMetal':7, 
                                  'metal':[-2.32, -1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(17.7828),
                                                       50),
                                  'metal_str':['m2.32', 'm1.71', 'm1.31', 'm0.71', 
                                               'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },
            'MILES_Padova00_kb_1.30_sub':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30_sub'),
                                  'nAges':25, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(15.8489),
                                                       25),
                                  'metal_str':['m1.71', 'm1.31', 'm0.71', 
                                               'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },  
            'MILES_Padova00_kb_1.30_sub_sub':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30_sub_sub'),
                                  'nAges':13, 'nMetal':5, 
                                  'metal':[-1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(15.8489),
                                                       13),
                                  'metal_str':['m1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },                   
           }

           
# ---| Some other useful bits & pieces |--------------------------------------------------     

# A dictionary to keep track of the filenames throughout the code. Default below will only
# be used to initalize the list on the first pass.
fn_list = {'snr_cube':None, # The rough measure of the SNR in the continuum and the emission lines
           'galdered_cube':None, # The raw datacubes corrected for Galactic extinction
           'lowess_pickle':None, # The pickle files generated during the LOWESS continuum fitting
           'lowess_cube':None, # The LOWESS fits to the continuum
           'ppxf_lams':None, # The file containing the various wavelength array (cropped, log-rebinned, etc ...) used by ppxf
           'ppxf_pickle':None, # The pickle files generated during the PPXF continuum fitting
           'ppxf_cube':None, # The fitted PPXF continuum
           'ppxf_delog_raw_cube':None, # The raw spectra after going through ppxf (log-rebin-delog-rebin). Useful to check how rebinning the ppxf fits affects the output.
           'ppxf_sol_map':None, # The PPXF "sol" array, i.e. the fitting parameters.
           'cont_mix_cube':None, # The cube containing the final mix of continuum fits (lowess vs ppxf)
           'elines_pickle':None, # The pickle files generated duirng the MPFIT fitting.
           'elines_spec_cube':None, # The emission line spectra
           'elines_params_cube':None, # The emission line parameters
           'elines_perror_cube':None, # The emission line parameter errors
           'elines_fit_status':None, # The emission line MPFIT fit status
           'ap_list':None, # The aperture list.
           'ap_map':None, # The aperture map.
           'ap_spec_cube':None, # The aperture "cube" 
           'dered_elines_params':None, # The extragalactic de-attenuated line parameters
           'dered_elines_perror':None, # The extragalactic de-attenuated line parameters errors
           'Av_map': None, # The extinction map.
           'fit_kin_pa':None, # A pickle file that contains the PA of the velocity map
           'ne_map':None, # The fits file containing the [SII] line ratio
           'pyqz_pickle': None, # The output of pyqz for each row, stored in a pickle file
          }

# ---| Custom colormaps for brutus |--------------------------------------------------------

# Here, I define the "alligator" colormap to give brutus a unique feel just because I can. 
# More importantly, I can decide to set the lower and upper bounds to B&W, so that spaxels
# outside the colorbar range are directly identifiable. Also make sure that nan's show up 
# in 50% grey. 

# Alligator: 6 nodes with linear evolution, suitable for B&W impression

cdict_alligator = {
'red'  :  ( (0.00, 0./255, 0./255),
            (0.00, 0./255, 0./255),      (0.2, 20./255., 20./255.),
            (0.40, 46./255., 46./255.),  (0.6, 108./255., 108./255.),
            (0.8, 207./255., 207./255.), (1.00, 255./255.,255./255.),
            (1.00, 255./255, 255./255),),

'green' : ( (0.00, 0./255, 0./255),
            (0.00, 25./255., 25./255.),  (0.2, 85./255., 85./255.),
            (0.4, 139./255., 139./255.), (0.6, 177./255., 177./255.),  
            (0.8, 234./255, 234./255),   (1.00, 248./255.,248./255.),
            (1.00, 255./255, 255./255),),
            
'blue':   ( (0.00, 0./255, 0./255),
            (0.00, 25./255., 25./255.),  (0.2, 81./255., 81./255.),
            (0.4, 87./255., 87./255.),   (0.6, 86./255., 86./255.),
            (0.8, 45./255, 45./255),     (1.00, 215./255.,215./255.),
            (1.00, 255./255, 255./255),),   
}

# Initiate the colorbar
alligator = plt.matplotlib.colors.LinearSegmentedColormap('alligator',cdict_alligator, 1024)
# Set the bad colors
alligator.set_bad(color=(0.5,0.5,0.5), alpha=1) 




      