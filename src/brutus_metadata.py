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

__version__ = '0.3.0'

# Where are we located ?
# Get the project directory from this file location !
brutus_dir = os.path.dirname(__file__)

# And get all the different useful locations associated
refdata_dir = os.path.join(brutus_dir,'../reference_data')

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


# ---| Information about the different stellar models for ppxf|---------------------------

sl_models ={'MILES_ppxf_default':{'sl_loc':os.path.join(refdata_dir,
                                                        os.path.join('stellar_templates',
                                                                     'miles_models')),
                                  'nAges':26, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(1),np.log10(17.7828),26),
                                  'metal_str':['m1.71','m1.31','m0.71','m0.40','p0.00',
                                               'p0.22'],
                                  'fwhm':2.51,
                                }, 

           }

           
# ---| Some other useful bits & pieces |--------------------------------------------------     

# A dictionary to keep track of the filenames throughout the code. Default below will only
# be used to initalize the list on the first pass.
fn_list = {'snr_cube':None, 
           'galdered_cube':None, 
           'lowess_pickle':None,
           'lowess_cube':None,
           'ppxf_lams':None,
           'ppxf_pickle':None,
           'ppxf_cube':None,
           'ppxf_delog_raw_cube':None,
           'ppxf_sol_map':None,
           'elines_pickle':None,
           'elines_spec_cube':None,
           'elines_params_cube':None,
           'elines_perror_cube':None,
           'elines_fit_status':None,
           'dered_elines_params':None,
           'dered_elines_perror':None,
           'ap_list':None,
           'ap_map':None,
           'ap_spec_cube':None,
           'Av_map': None,
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




      