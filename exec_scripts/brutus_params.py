# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains some global variables and other metadata used throughout the
 brutus set of Python modules to fit IFU datacubes. 

 This is where (if all goes according to plan) users can modify the different
 fit parameters to suit their needs. 

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import pickle
import os
from brutus_metadata import * # Load some ref. emission lines

# What si the name of the resulting pickle file going to be ?
pickle_fn = 'HCG91c_params.pkl'

# Where is everything ?
# Get the project directory from the file location !
proj_dir = os.path.dirname(__file__)

brutus_params = {

    # ---| General variables |------------------------------------------------------------
    'target': 'HCG91c', # The generic name for files created by brutus
    'z_target': 0.023983, # Guess redshift of galaxy
    'inst': 'MUSE', # Which instrument took the data ?
    'multiprocessing': True,
    'verbose':True,
    'warnings':'ignore', # 'ignore' = silence known (?) warnings about all-nan's spaxels
                         # 'default' = default system behavior.

	# ---| Location and name of the data file to process |--------------------------------
    'data_loc': os.path.join(proj_dir,'reduc/'), # relative path from this file loc !
    'data_fn': 'WFM-HCG91c-3-0_DATACUBE_FINAL_WFM-NOAO-N_475-935.fits',

	# ---| Generic locations for plots and products |-------------------------------------
    'plot_loc': os.path.join(proj_dir,'plots/') ,
    'prod_loc': os.path.join(proj_dir,'products/'),
    'tmp_loc': os.path.join(proj_dir,'products/tmp/'),
    
    'fn_list_fn': 'filenames_dictionary.pkl', # Name of the dictionary for filenames
    
    # ---| Constructing SNR maps |--------------------------------------------------------
    'cont_range':[6400.,6500.], # The range over which to get the continuum SNR.
    
    # ---| Correcting for the Galactic extinction |---------------------------------------
    'gal_curve': 'f99', # Set to 'f99' to follow NED.
    'gal_rv': 3.1, # Set to 3.1 to follow NED.
    'Ab': 0.069, # Ab (Landolt) from NED.
    'Av': 0.052, # Av (Landolt) from NED.
    
    # ---| Continuum fitting |------------------------------------------------------------
    # Lowess specific variables
    'lowess_snr_min': 0,  # What's the lowest SNR to run the LOWESS fit ?
    'lowess_snr_max': None, # What's the highest SNR to run the LOWESS fit ? None = max
    'lowess_it': 5, # Number of iteration for sigma-clipping to get rid of outliers
    'lowess_frac': 0.05, # % of the array used for deriving each point. 0.05 = sweet spot?
    
    # PPXF specific variables
    'ppxf_snr_min': 3,  # What's the lowest SNR to run the PPXF fit ?
    'ppxf_snr_max': None, # What's the highest SNR to run the PPXF fit ? None = max
    'ppxf_sl_name': 'MILES_ppxf_default',   # Name of the stellar libraries for ppxf
    'ppxf_sampling':1, # Sampling of the spectra when log-rebining them.
    # ppxf_sampling > 1 can help reduce errors when log-rebining and de-log-rebining the 
    # ppxf fit output. But ppxf_sampling>1 is very expansive time-wise. Worth it only with 
    # very good quality continuum - and/or and powerful computer.
    'ppxf_moments':2, # Moments for the line profile (2 = v,sigma, 4 = v, sigma, h3, h4)
    'ppxf_degree':-1, # Degree of additive polynomials (-1 = None)
    'ppxf_mdegree':10, # Degree of multiplicative polynomials
    'ppxf_regul_err':0.004, # Desired regularization error
    
    # ---| Emission line fitting | -------------------------------------------------------
    # How do I want to subtract the continuum: give the range of snr (in the continuum!)
    # and the corresponding technique.
    'which_cont_sub':{'0->max':'lowess'},
    
    # For the inital delta-v guess, the code looks at one single line. Which one is the
    # strongest in your case ? Give here the un-redshifted wavelength.
    'ref_dv_line':ha,
    
    # What kind of profile for emission lines ?
    # 'gauss_hist' = gaussian, accounting for bin size (the correct way).
    'line_profile':'gauss', # 'gauss', 'gauss-herm'
    
    'elines_snr_min': 0,  # What's the lowest SNR to run the line fit ?
    'elines_snr_max': None, # What's the highest SNR to run the line fit ? None = max
    
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
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':''}, # Velocity
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':''}, # Dispersion >0
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
                  ],
              'b':[[n2h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'},
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
                  ],
              'c':[[n2l,0,20],
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./3.*p(b)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'},
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
                  ],
              'd':[[s2h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'}, 
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
                  ],
              'e':[[s2l,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'}, 
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
                  ],
              'f':[[hb,0,20],
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':''},
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':''},
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
                  ],
              'g':[[o3h,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(f)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(f)'},
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
                  ],
              'h':[[o3l,0,20], 
                   {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./2.98*p(g)'},
                   {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(f)'}, 
                   {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(f)'},
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
                   {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
                  ],
             },
             
    # ---| Structure detection and aperture extraction |----------------------------------
    # What line flux do I want to use to extract my aperture map ? Use a list of the 
    # elines keys for identification.
    'ap_map_lines':['a','g'],
    
    # The default radius for the apertures in the automated peak detection routine, in pixel.
    # Choosing a value corresponding to the seeing is a good first guess.
    'ap_radius':3.0,             
    
    # ---| Extragalactic reddening/attenuation correction |------------------------------- 
    'hahb_0': 2.85, # theoretical reference flux ratio of Halpha/Hbeta. 2.85 for Case B 
                    # recombination, 3.1 for AGNs.
    'egal_curve': 'fd05', # Reddening/attenuation law to use. Valid values are:
                          # 'fd05': Fischera & Dopita (2005) use fv=3.08 and rva=4.5 for a 
                          #         curve similar to Calzetti (2000).
                          # 'cal00': Calzetti (2000), use with rv=4.05 for default.
                          # 'cal00*': Calzetti (2000), with 0.44 correction factor.
                          # 'f99': Fitzpatrick (1999)
                          # 'ccm89': Cardelli, Clayton & Mathis (1989)
    'egal_rv': 3.08, # The value of Rv. Use 3.08 for a Calzetti-like law, when using 'fd05'.
    'egal_rva':4.3,  # The value of Rva, if using 'fd05'.       
    }
    

# Export this as a pickle file    
f1 = open(pickle_fn, 'w')
pickle.dump(brutus_params, f1)
f1.close()            