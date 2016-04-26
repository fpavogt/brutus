# -*- coding: utf-8 -*-
#
# This file contains some global metadata used throughout the
# BRIAN code to fit MUSE data cubes. 
# 
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
import os

# ---| Some basic parameters |------------------------------------------------------------

__version__ = '0.1.1'

# Where are we located ?
# Get the project directory from this file location !
brian_dir = os.path.dirname(__file__)

# And get all the different useful locations associated
refdata_dir = os.path.join(brian_dir,'reference_data')

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


# ---| Information about the different stellar models |-----------------------------------
'''
sl_models ={'MILES_ppxf_default':{'sl_loc':os.path.join(ppxf_dir,'miles_models'),
                                  'nAges':26, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(1),np.log10(17.7828),26),
                                  'metal_str':['m1.71','m1.31','m0.71','m0.40',
                                               'p0.00','p0.22'],
                                  'fwhm':2.51,
                                }, 

           }
'''
           
# ---| Some other useful bits & pieces |--------------------------------------------------     

# A dictionary to keep track of the filenames throughout the code. Default below will only
# be used to initalize the list on the first pass.
fn_list = {'snr_cube':None,
           'lowess_pickle':None,
           'lowess_cube':None,
           'ppxf_pickle':None,
           'ppxf_cube':None,
           'elines_pickle':None,
           'elines_spec_cube':None,
           'elines_params_cube':None,
           'elines_perror_cube':None,
           'elines_fit_status':None,
           'ap_list':None,
           'ap_map':None,
           'ap_spec_cube':None,
          }
      