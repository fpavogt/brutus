# -*- coding: utf-8 -*-
#
# This program is designed to fit some emission lines in a given MUSE datacube, after
# removing the continuum, either via a polynomial fitting (noisy case) or using ppxf
# (if the SNR allows to do so semi-meaningfully).
#
# This code is specifically designed to handle: 
#   a) the varying spectral resolution of MUSE data as a function of wavelength, and 
#   b) the residual error in the emission line centroid due to difficult LSF issues in the
#      data reduction procedure.
#
# The idea for point a) is to properly convolve the template spectra, and for b) to fit 
# the emission lines with a global reference velocity plus a small (<~5km/s) individual 
# velocity offset. 
#
#
# TO DO: 
# - understand the output of ppxf, i.e. be able to remove the continuum from the raw
#   spectra, before doing the emission line fitting.
# - propagate variance from MUSE to ppxf (currently just using a constant number)
# - implement emission line fitting routine for MUSE.
#
# - implement a source detection+extraction routine for getting the integrated light from
#   HII regions.
#
# Created February 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# ---------------------------| Import the packages |--------------------------------------
import pickle
import sys
import brian

import datetime

# ----------------------------------------------------------------------------------------
#---------------------------| Set up the scene, get the params, etc ... |-----------------

start_time = datetime.datetime.now()

# Get name of metadata file - fed to the code when starting
try:
    params_fn = sys.argv[1]
    f1 = open(params_fn, 'r')
    bp = pickle.load(f1)
    f1.close()
except:
    sys.exit('Failed to load pickle file with parameters. Have you supplied any ?')


# ----------------------------------------------------------------------------------------
# --------------------------| The different processing steps to be run |------------------
proc_steps = [
    # -----
    # First, compute some SNR maps for the continuum and Halpha ...
    {'step':'snr_maps', 'run':False, 'suffix':'00', 
     'args':{'do_plot':True, # Do I want to save some pretty pictures ?
            }
    },
    # -----
    # First, fit the continuum using a LOWESS filter. Much more elegant than a polynomial!
    {'step':'fit_continuum', 'run':False, 'suffix':'01',
     'args':{'start_row': 9, # Where to start the fitting ? None = 0
             'end_row': 9, # Where to end the fitting ? None = max
             'method': 'lowess', # What kind of fitting is this ?
             },
    },
    # Then, we do it again using PPXF, but only for the spaxels with "decent" SNR ...
    #{'step':'fit_continuum', 'run':False, 'suffix':'02',
    # 'args':{'start_row': 150, # Where to start the fitting ? None = 0
    #         'end_row': 150, # Where to end the fitting ? None = max
    #         'method': 'ppxf', # What kind of fitting is this ?
    #         }
    #},
    # -----
    # Let's move on to emission line fitting.
    {'step':'fit_elines', 'run':False, 'suffix':'03',
     'args':{'start_row':0, # Where to start the fitting ? None = 0
             'end_row':None, # Where to end the fitting ? None = max
             },
    },
    # Construct a datacube with the emission line parameters, using the output of the 
    # previous step
    {'step':'make_elines_cube', 'run':False, 'suffix':'04',
     'args':{'prev_suffix':'03'},
    },
    {'step':'plot_elines_cube', 'run':True, 'suffix':'05',
     'args':{'prev_suffix':'04', 'vrange':[7220,7400], 'sigrange':[45,60]},
    },
    # Now, create correct for the reddening of the spectra
    #{'step':'dered', 'run':False, 'suffix':'05',
    # 'args':{},
    #},
    # Compute the electron density
    #{'step':'electron_density', 'run':False, 'suffix':'05',
    # 'args':{},
    #},
    # Derive the oxygen abundance and ionization parameters using pyqz (separate install)
    #{'step':'get_QZ', 'run':False, 'suffix':'05',
    # 'args':{},
    #},
    
    ]
    

# ----------------------------------------------------------------------------------------
# -------------| Call the individual master step functions as required |------------------

prev_suffix = None
for step in proc_steps:
    step_name   = step['step']
    step_run    = step['run']
    step_suffix = step['suffix']
    step_args   = step['args']
    func_name = 'run_'+step_name
    func = getattr(brian,func_name)
    if step_run:
        func(bp,
             suffix = step_suffix,
             **step_args)
    if step_suffix != None:
        prev_suffix = step_suffix

 
# ------------------------| All done ! |--------------------------------------------------
duration = datetime.datetime.now() - start_time
print 'All done in %.01f seconds.' % duration.total_seconds() 
 
# -----------------------| End of the World as we know it |-------------------------------