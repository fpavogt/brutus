# -*- coding: utf-8 -*-
#
# This file contains some global variables and other metadata used throughout the
# BRIAN (BRIlliant AcroNym) routine to fit MUSE spectra. 
#
# This is where (if all goes according to plan) users can modify the different
# fit parameters to suit their needs. 
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import pickle
import os
from brian_metadata import * # Load some ref. emission lines


# What do I want to name the resulting pickle file ?
pickle_fn = 'HCG91c_params.pkl'

# Where is everything ?
# Get the project directory from the file location !
proj_dir = os.path.dirname(__file__)

brian_params = {

    # ---| General variables |------------------------------------------------------------
    'target': 'HCG91c', # The generic name for files created by BRIAN
    'z_target': 0.023983, # Guess redshift of galaxy
    'inst': None, # Which instrument took the data ?
    'multiprocessing': True,
    'verbose':True,

	# ---| Location and name of the data file to process |--------------------------------
    'data_loc': os.path.join(proj_dir,'../reduc/'), # relative path from this file loc !
    'data_fn': 'WFM-HCG91c-3-0_DATACUBE_FINAL_WFM-NOAO-N_475-935.fits',

	# ---| Generic locations for plots and products |-------------------------------------
    'plot_loc': os.path.join(proj_dir,'plots/') ,
    'prod_loc': os.path.join(proj_dir,'products/'),
    'tmp_loc': os.path.join(proj_dir,'products/tmp/'),
    
    # ---| Constructing SNR maps | -------------------------------------------------------
    'cont_range':[6400.,6500.], # The range over which to get the continuum SNR.
    
    # ---| Continuum fitting | -----------------------------------------------------------
    # Lowess specific variables
    'lowess_snr_min': 0,  # What's the lowest SNR to run the LOWESS fit ?
    'lowess_snr_max': None, # What's the highest SNR to run the LOWESS fit ? None = max
    
    # PPXF specific variables
    'ppxf_snr_min': 20,  # What's the lowest SNR to run the PPXF fit ?
    'ppxf_snr_max': None, # What's the highest SNR to run the PPXF fit ? None = max
    'sl_name': 'MILES_ppxf_default',   # Name of the stellar libraries for ppxf
    
    
    # ---| Emission line fitting | -------------------------------------------------------
    # How do I want to subtract the continuum: give the range of snr (in the continuum!)
    # and the corresponding technique.
    'which_cont_sub':{'0->max':'lowess'},
    
    # For the inital delta-v guess, the code looks at one single line. Which one is the
    # strongest in your case ? Give here the un-redshifted wavelength.
    'ref_dv_line':ha,
    
    # What kind of profile for emission lines ?
    # 'gauss_hist' = gaussian, accounting for bin size (the correct way).
    'line_profile':'gauss_hist',
    
    # What lines do I want to fit, and how ?
    # This is very close to defining the mpfit "parinfo" structure, but not quite.
    # For each line, first provide the reference wavelength, the offset from the "main"
    # component (NOT the total offset - used when fitting multiple gaussians to 1 line),
    # and the guess velocity dispersion (in km/s).
    # Then, individual dictionaries for the intensity, velocity and sigma of each lines
    # allow to define "fixed" parameters (that will NOT befitted), set lower/upper bounds,
    # and tie a parameters with another one (e.g. line intensity, velocity, etc ...)
    
    'elines':{'a':[[ha,0,20], # Line reference wavelength, guess dispersion (km/s)
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, # Intensity >0
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # Velocity
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':''}, # Dispersion >0
                  ],
              'b':[[n2h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(a)'}, 
                  ],
              'c':[[n2l,0,20],
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./3.*p(b)'}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(a)'}, 
                  ],
              'd':[[s2h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(a)'}, 
                  ],
              'e':[[s2l,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(a)'}, 
                  ],
              'f':[[hb,0,20],
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':''},
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':''}, 
                  ],
              'g':[[o3h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(f)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(f)'}, 
                  ],
              'h':[[o3l,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./2.98*p(g)'},
                   {'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':'p(f)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[0.,100.], 'tied':'p(f)'}, 
                  ],
             },
    }
    

# Export this as a pickle file    
f1 = open(pickle_fn, 'w')
pickle.dump(brian_params, f1)
f1.close()            