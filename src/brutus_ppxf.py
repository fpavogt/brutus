# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains functions to call ppxf to fit the continuum in a given spectra.
 It is based on the examples provided by M. Cappellari within ppxf istelf
 (e.g. ppxf_kinematics_example_sauron.py & ppxf_population_example_sdss.py), but has been
 modified to account for a wavelength dependant spectral resolution of the data. 

 Basically, this is designed with MUSE in mind.

 Created February 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

#from __future__ import print_function

from astropy.io import fits as pyfits
from scipy import ndimage
import numpy as np
import glob
import matplotlib.pyplot as plt
from time import clock
import datetime
import warnings

import os
import scipy.signal as signal
import sys

from ppxf import ppxf
import ppxf_util as util

import brutus_tools
from brutus_metadata import *
        
# ----------------------------------------------------------------------------------------
def sl_resample(fn, z=0, inst = 'MUSE', sl_name='MILES_ppxf_default', do_cap=True):
    '''
    This function resamples a given spectrum (i.e. from a stellar library)
    to match theresolution of the instrument that took the data.
    
    :param: fn: filename of stellar template
    :param: z: redshift of the target observations
    :param: inst: which instrument are we using
    :param: sl_name: which templates do we use ?
    :do_cap: do the convolution using M.Cappellari in-built function
    
    Returns the convolved spectrum
    '''
        
    if not(os.path.isfile(fn)):
        sys.exit('"%s": not a valid file.' % fn)
    
    # Load the ssp spectra
    hdu = pyfits.open(fn)
    ssp_raw = hdu[0].data
    ssp_header = hdu[0].header
    ssp_out = np.zeros_like(ssp_raw)
    hdu.close()
    
    # The wavelength step and wavelength array
    dlam = ssp_header['CDELT1']
    lams = np.arange(0,ssp_header['NAXIS1'],1)*ssp_header['CDELT1'] + \
                      ssp_header['CRVAL1']
    
    # What is the target resolution I want to achieve ?
    if inst:
        R_target = brutus_tools.inst_resolution(inst=inst, get_ff=False, show_plot=False)
        fwhm_target = lams*(z+1)/R_target(lams*(z+1))
    else:
        # Assuming MUSE - probably a bad idea that I will forget about ... Issue a warning
        # to be on the safe side.
        warnings.warn('No instrument specified. Assuming MUSE for the ppxf input spectra')
        R_target = brutus_tools.inst_resolution(inst='MUSE', get_ff=False, show_plot=False)
        fwhm_target = lams*(z+1)/R_target(lams*(z+1))

    # WARNING: what if the resolution of the models is NOT sufficient ?
    if np.any((fwhm_target - sl_models[sl_name]['fwhm']) < 0):
        sys.exit('WARNING: stellar templates have too high fwhm given the redshift '+
                      'of the target ! Ignoring z for now ... you really should use '+
                      'other models !') 
    
	# Use the ppxf-in-built gaussian_fliter1d home made by M. Cappellari (10%faster !)
    if do_cap:
       
       # Compute the correction in FWHM to be applied.	        
       fwhm_diff = np.sqrt(fwhm_target**2 - sl_models[sl_name]['fwhm']**2) 
       # Transform this from Angstroem to resolution elements
       fwhm_diff = fwhm_diff/dlam # In resolution elements 
       # Compute the corresponding sigma for the gaussian profile                      
       sigma = (2*np.sqrt(2*np.log(2)))**-1 * fwhm_diff
       
       # Feed this to M. Cappellari elegant routine. N.A.: I checked that function, and 
       # I can reproduce its results (see below) doing things my way.    
       ssp_out = util.gaussian_filter1d(ssp_raw,sigma) 
       
    else: # do things my way - this is perfectly equivalent, but slower ...
	        
       # What is the max FWHM in resolution elements ?
       fwhm_pix_max = np.max(fwhm_target/dlam)
        
       # Very well, let's use this to define my box size. Make it 8 times
       # the max FWHM to correctly include the wings of the gaussian kernel.
       window_size = 8*np.ceil(fwhm_pix_max)
       
       # Make sure the window size is odd (does it matter ?)
       if window_size % 2 ==0:
          window_size +=1
       
       window_size = np.int(window_size)
       
       # Ok, start the convolution. I do it semi-manually, for cases where it
       # varies as a function of wavelength (i.e. MUSE)
       for i in range(len(ssp_raw)):     
          this_lam = lams[i] # The wavelength of the i-th element
          this_fwhm_target = fwhm_target[i]
           
          # Compute the correction in FWHM to be applied. 
          fwhm_diff = np.sqrt(this_fwhm_target**2 - sl_models[sl_name]['fwhm']**2)
          # Transform this from Angstroem to resolution elements  
          fwhm_diff = fwhm_diff/dlam 
          # Compute the corresponding sigma for the gaussian profile                   
          sigma = (2*np.sqrt(2*np.log(2)))**-1 * fwhm_diff # In res. elements
          
          # Here's my gaussian kernel   
          w = signal.gaussian(window_size,sigma,sym = True)
        
          # If the window lands outside the spectrum, do nothing ...      
          if (i < window_size/2.) or (i > len(ssp_raw) - window_size/2.): 
              ssp_out[i] = 0.0        
          else:
              # Here, I am using correlate, because it is exactly what
              # np.convolve does. It's a way to get the same numbers out,
              # (down to the machine precision limit).
              ssp_out[i] = np.correlate(w/w.sum(),
                                    ssp_raw[i-(window_size-1)/2:
                                                  i+(window_size-1)/2+1],
                                    mode='valid')               
              '''
              # Some comparison with similar methods, to see if I 
              # really understand what I do ...
                 
              out_spec2 = np.convolve(w/w.sum(),
                                  	raw_spec[:,1][i-(window_size-1)/2:
                                                  i+(window_size-1)/2+1],
                                    mode='valid')
              # -> identical. Makes sense, since np.convolve wraps around
              # np.correlate.
                                        
              # Next, check the "official" gaussian_filter1d                      
              out_spec3 = ndimage.gaussian_filter1d(raw_spec[:,1][i-(window_size-1)/2:
                                                                i+(window_size-1)/2+1],
                                                      sigma)
              # -> different at 10**-6 level. Results are different 
              # because I have a larger window size compared to what 
              # gaussian_filter1d does by default.                                     
                                                      
              # But checking closer with the function inside
              # gaussian_filter1d, the results are actually identical, 
              # when I use the same weights !
              out_spec4 = ndimage.filters.correlate1d(raw_spec[:,1][i-(window_size-1)/2:
                                                    i+(window_size-1)/2+1],
                                                    w/w.sum())                                     
              import pdb
              pdb.set_trace()
                 
              # Bottom line I'm doing things correctly. Slightly more 
              # accurately than in the example of M. Cappellari.
              '''                                                                                          
    '''                                    
    plt.close(1)
    plt.figure(1)      
    plt.plot(lams,ssp_raw,'k-')
    plt.plot(lams,ssp_out,'r-')
    plt.show()
    '''   
    return ssp_out

def setup_spectral_library(velscale, inst, sl_name, fit_lims, z):
    '''
    This prepares the set of stellar templates to be fed to ppxf.
    
    :param: velscale - the velocity scale for the log-rebinning
    :param: inst - the instrument that took the data, e.g. 'MUSE'
    :param: sl_name - the name of the stellar library to use
    :param: fit_lims - the spectral range for the fit (data & templates overlap)
    :param: z - redshift of target - needed for correct reshaping of fwhm_target.
    
    :return: ...lots of stuf...
    '''
    
    # Fetch the stellar template files
    fns = glob.glob(os.path.join(sl_models[sl_name]['sl_loc'],'*'))
    fns.sort()
    
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the data spectra, to determine
    # the size needed for the array which will contain the template spectra.
    hdu = pyfits.open(fns[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lams_temp = np.arange(0,h2['NAXIS1'],1)*h2['CDELT1'] + h2['CRVAL1']
    
    mask = (lams_temp >= fit_lims[0]) & (lams_temp <= fit_lims[1])
    ssp = ssp[mask]
    lams_temp = lams_temp[mask]
    lam_range_temp = np.array([lams_temp[0],lams_temp[-1]])
    
    log_lam, ssp_new, velscale = brutus_tools.log_rebin(lam_range_temp, ssp, 
                                                        velscale=velscale)
    
    # Create a three dimensional array to store the two dimensional grid of model spectra
    nAges = sl_models[sl_name]['nAges']
    nMetal = sl_models[sl_name]['nMetal']
    templates = np.empty((ssp_new.size, nAges, nMetal))

    # These are the array where we want to store the characteristics of each SSP model
    logAge_grid = np.empty((nAges, nMetal))
    metal_grid = np.empty((nAges, nMetal))

    # These are the characteristics of the adopted rectangular grid of SSP models
    logAge = sl_models[sl_name]['logAge']
    metal = sl_models[sl_name]['metal']

    # Here we make sure the spectra are sorted in both [M/H]
    # and Age along the two axes of the rectangular grid of templates.
    metal_str = sl_models[sl_name]['metal_str']
    for k, mh in enumerate(metal_str):
        files = [s for s in fns if mh in s]
        for j, filename in enumerate(files):
            # Resample the spectra to match the instrumental resolution
            ssp = sl_resample(filename, z=z, inst=inst, sl_name=sl_name, do_cap=True)
            ssp = ssp[mask]
            # Log-rebin the spectrum
            log_lam, ssp_new, velscale = brutus_tools.log_rebin(lam_range_temp, ssp, 
                                                                velscale=velscale)
            
            templates[:, j, k] = ssp_new  # Templates are *not* normalized here
            logAge_grid[j, k] = logAge[j]
            metal_grid[j, k] = metal[k]

    return templates, lam_range_temp, logAge_grid, metal_grid, log_lam

#------------------------------------------------------------------------------

def ppxf_population(specerr, templates=None, velscale=None, start = None, 
                    goodpixels = None, plot=False, moments=4, degree=-1, vsyst=None, 
                    clean=False, mdegree=10, regul=None):
    '''
    Function that calls ppxf, designed to be called from multiprocessing (i.e. only 
    require 1 argument). Everything else is taken care of outside of this.
    '''
    
    galaxy = specerr[0]
    noise = specerr[1]
    
    # Only run ppxf if I have a non-nan's spectrum. Even one nan's breaks ppxf ... sigh.
    if not(np.any(np.isnan(galaxy))):

        pp = ppxf(templates, galaxy, noise, velscale, start,
                  goodpixels=goodpixels, plot=False, moments=moments, degree=degree,
                  vsyst=vsyst, clean=clean, mdegree=mdegree, regul=regul, quiet=True)
        
        # In a  perfect world, I would return pp. But that is a very massive variable.
        # 1 row = 300 spaxels = 5GB => 1 MUSE cube = 1.5TB !!!
        # Instead, just store what I need/trust. The best fit, the galaxy spectra (for 
        # sanity checks) and the sol array (with the kinematics information, that I 
        # sort-of trust.    
        return [pp.bestfit, pp.galaxy, pp.sol]
    
    else:
        return None
    

#------------------------------------------------------------------------------
