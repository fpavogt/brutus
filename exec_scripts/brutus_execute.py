# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This program is designed to fit some emission lines in a given IFU datacube, after
 removing the continuum, either via a polynomial fitting (noisy case) or using ppxf
 (if the SNR allows to do so semi-meaningfully).

 This code is specifically designed to handle: 
   a) the varying spectral resolution of MUSE data as a function of wavelength, and 
   b) the residual error in the emission line centroid due to difficult LSF issues in the
      data reduction procedure.

 The idea for point a) is to properly convolve the template spectra, and for b) to fit 
 the emission lines with a global reference velocity plus a small (<~5km/s) individual 
 velocity offset. 


 TO DO: 
 - propagate error when correcting for extragalctic reddening

 Created February 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# ---------------------------| Import the packages |--------------------------------------
import pickle
import sys
import os
import brutus
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

# Is there a dictionary of filenames already in place ? If not, create one
fn =os.path.join(bp['prod_loc'],bp['fn_list_fn'])
if not(os.path.isfile(fn)):
    if bp['verbose']:
        print '-> This looks like a fresh run. Creating the dictionary of filenames.'
    from brutus_metadata import fn_list
    f = open(fn,'w')
    pickle.dump(fn_list,f)
    f.close()

# ----------------------------------------------------------------------------------------
# --------------------------| The different processing steps to be run |------------------
proc_steps = [
    # ----- Setting up the scene ---------------------------------------------------------
    # First, compute some SNR maps for the continuum and emission lines ...
    {'step':'snr_maps', 'run':False, 'suffix':'00', 
     'args':{'do_plot':True, # Do I want to save some pretty pictures ?
            },
    },
    # Correct for galactic reddening via the values given on NED
    {'step':'gal_dered', 'run':False, 'suffix':'01',
     'args':{'do_plot':True},
    },
    # ----- Continuum fitting ------------------------------------------------------------
    # First, fit the continuum using a LOWESS filter. Much more elegant than a polynomial!
    {'step':'fit_continuum', 'run':False, 'suffix':'02',
     'args':{'start_row': 0, # Where to start the fitting ? None = 0
             'end_row': 0, # Where to end the fitting ? None = max
             'method': 'lowess', # What kind of fitting is this ?
            },
    },
    # Gather the output in a single datacube
    {'step':'make_continuum_cube', 'run':False, 'suffix':'03',
     'args':{'method':'lowess',
            },
    },
    # Then, we do it again using PPXF, but only for the spaxels with "decent" SNR ...
    # (ppxf = separate install)
    {'step':'fit_continuum', 'run':False, 'suffix':'04',
     'args':{'start_row': 0, # Where to start the fitting ? None = 0
             'end_row': 0, # Where to end the fitting ? None = max
             'method': 'ppxf', # What kind of fitting is this ?
             }
    },
    # Gather the output in a single datacube
    {'step':'make_continuum_cube', 'run':False, 'suffix':'05',
     'args':{'method':'ppxf',
            },
    },
    # Make some plots of the output of ppxf
    {'step':'plot_ppxf_sol', 'run':False, 'suffix':'06',
     'args':{'vrange':[7210,7400], # Velocity range in km/s for the colorbar 
             'sigrange':[20,80],
            },
    },
    # ----- Emission line fitting --------------------------------------------------------
    # Let's move on to emission line fitting.
    {'step':'fit_elines', 'run':False, 'suffix':'07',
     'args':{'start_row':0, # Where to start the fitting ? None = 0
             'end_row':0, # Where to end the fitting ? None = max
            },
    },
    # Construct a datacube with the emission line parameters, using the output of the 
    # previous step
    {'step':'make_elines_cube', 'run':False, 'suffix':'08',
     'args':{},
    },
    {'step':'plot_elines_cube', 'run':False, 'suffix':'09',
     'args':{'vrange':[7210,7400], # Velocity range in km/s for the colorbar 
             'sigrange':[40,60], # Dispersion range in m/s for the colorbar
            }, 
    },
    # ----- Fit inspection ---------------------------------------------------------------
    # Interactive window displaying the fit, and comparison between the different 
    # continuum removal.
    {'step':'inspect_fit','run':False, 'suffix':'10',
     'args':{'irange':[1,1e4], 
             'vrange':[7200,7400],
            },
    },
    # ----- Structures and apertures business --------------------------------------------
    # Detect structures automatically and possibly refine the selection manually
    {'step':'find_structures', 'run':False, 'suffix':'11',
     'args':{'automatic_mode':True, # Detect structures automatically ?
             'interactive_mode':True, # Refine structure selection manually ?
            },
    },
    # Extract the integrated spectra in the apertures, and save as a new cube
    {'step':'make_ap_cube', 'run':False, 'suffix':'12',
     'args':{'do_plot':True,
            },
    },
    # Now, you probably want to re-process this aperture file, right ? 
    # I suggest to not add extra steps below, but rather create a dedicated copy of
    # brutus_params.py, and treat the aperture cube as a new cube.
    # This approach is much less invasive, and somewhat more elegant. After all, the 
    # aperture cube is build to be processed on a spaxel-by-spaxel basis ! It's somewhat
    # time consuming, but you'll have to live with it for now.
    # ----- Emission line analysis -------------------------------------------------------
    # Now, correct for the extragalactic reddening of the spectra, using Ha and Hb
    {'step':'extragal_dered', 'run':True, 'suffix':'13',
     'args':{'do_plot':True},
    },
    # Compute the electron density
    #{'step':'electron_density', 'run':False, 'suffix':'05',
    # 'args':{},
    #},
    # Derive the oxygen abundance and ionization parameters using pyqz 
    # (pyqz = separate install)
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
    # Call a function based on a string ... pretty sweet !
    func = getattr(brutus,func_name) 
    if step_run:
        # Here, I want to maintain a dictionary of filenames, to be used accross functions
        # For each step, load the dictionary, feed it to the function, and save it update
        # it when it's ll done. Each unction returns the updated dictionary !
        fdf = os.path.join(bp['prod_loc'],bp['fn_list_fn'])
        f = open(fdf,'r')
        fn_list = pickle.load(f)
        f.close()
        
        # Launch the function
        fn_list = func(fn_list, bp, suffix = step_suffix, **step_args)
        
        # Save the updated dictionary of filenames
        f = open(fdf,'w')
        fn_list = pickle.dump(fn_list,f)
        f.close()

# ------------------------| All done ! |--------------------------------------------------
duration = datetime.datetime.now() - start_time
print 'All done in %.01f seconds.' % duration.total_seconds() 
 
# -----------------------| End of the World as we know it |-------------------------------