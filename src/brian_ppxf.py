# -*- coding: utf-8 -*-
#
# This routine calls ppxf to fit the continuum in a given spectra.
# It is based on the examples provided by M. Cappellari within ppxf istelf
# (e.g. ppxf_kinematics_example_sauron.py & ppxf_population_example_sdss.py), but has been
# modified to account a wavelength dependant spectral resolution of the data. 
#
# Basically, this is designed with MUSE in mind.
#
# Created February 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

from __future__ import print_function

try:
    import pyfits
except ImportError:
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

import elf_tools
from elf_metadata import *
        
# ----------------------------------------------------------------------------------------
def sl_resample(fn, z=0, inst = 'MUSE', sl_name='MILES_ppxf_default', do_cap=True):
    '''
    This function resamples a given spectrum (i.e. from a stellar library)
    to match the MUSE resolution varying with wavelength.
    
    :param: fn: filename of stellar template
    :param: z: redshift of MUSE observations
    :param: inst: which instrument are we using
    :param: sl_name: which templates do we use ?
    :do_cap: do the convolution using M.Cappellari in-built function
    
    Returns the convolved spectrum
    '''
    
    if not(inst in ['MUSE']):
        sys.exit('Instrument unsupported ...')
    
    if not(sl_name in sl_models.keys()):
        sys.exit('Stellar library not supported ...')
        
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
    R_target = elf_tools.inst_resolution(inst=inst, get_ff=False, show_plot=False)
    fwhm_target = lams/R_target(lams)/(z+1)

    # WARNING: what if the resolution of the models is NOT sufficient ?
    if np.any((fwhm_target - sl_models[sl_name]['fwhm']) < 0):
        warnings.warn('WARNING: stellar templates have too high fwhm given the redshift '+
                      'of the target ! Ignoring z for now ... you really should use '+
                      'other models !') 
        sys.stdout.flush()                  
        fwhm_target = lams/R_target(lams)#/(z+1) 
    
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
       
    else: # do things my way - this is perfectly equivalent, but slower ?
	        
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
    
    import pdb
    pdb.set_trace()
    '''   
    return ssp_out

def setup_spectral_library(velscale, inst, sl_name, fit_lims, z):
    '''
    This prepares the set of stellar templates to be fed to ppxf.
    
    :param: velscale: the velocity scale for the log-rebinning
    :param: fwhm_gal: the fwhm of the target galaxy spectra - a poly1d structure !
    :param: sl_loc: the location of the stellar templates to be used.
    :param: fit_lims: the spectral range for the fit.
    :param: z: redshift of target - needed for correct reshaping of fwhm_target.
    
    :return: ...lots of stuf...
    '''
    
    fns = glob.glob(os.path.join(sl_models[sl_name]['sl_loc']+'*'))
    fns.sort()
    
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SDSS galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    #
    hdu = pyfits.open(fns[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lams_temp = np.arange(0,h2['NAXIS1'],1)*h2['CDELT1'] + h2['CRVAL1']
    
    mask = (lams_temp >= fit_lims[0]) & (lams_temp <= fit_lims[1])	
    ssp = ssp[mask]
    lams_temp = lams_temp[mask]
    lamRange_temp = np.array([lams_temp[0],lams_temp[-1]])
    
    sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, velscale=velscale)
	
    # Create a three dimensional array to store the
    # two dimensional grid of model spectra
    #
    nAges = sl_models[sl_name]['nAges']
    nMetal = sl_models[sl_name]['nMetal']
    templates = np.empty((sspNew.size, nAges, nMetal))

    # These are the array where we want to store
    # the characteristics of each SSP model
    #
    logAge_grid = np.empty((nAges, nMetal))
    metal_grid = np.empty((nAges, nMetal))

    # These are the characteristics of the adopted rectangular grid of SSP models
    #
    logAge = sl_models[sl_name]['logAge']
    metal = sl_models[sl_name]['metal']

    # Here we make sure the spectra are sorted in both [M/H]
    # and Age along the two axes of the rectangular grid of templates.

    metal_str = sl_models[sl_name]['metal_str']
    for k, mh in enumerate(metal_str):
        files = [s for s in fns if mh in s]
        for j, filename in enumerate(files):
            # Resample the spectra
            #start_time = datetime.datetime.now()
            
            ssp = sl_resample(filename, z=z, inst=inst, sl_name=sl_name, do_cap=True)
            
            #dt =datetime.datetime.now() - start_time
            #print (dt.total_seconds())
            
            ssp = ssp[mask]
            
            sspNew, logLam2, velscale = util.log_rebin(lamRange_temp, ssp, 
                                                       velscale=velscale)
                                                                 
            templates[:, j, k] = sspNew  # Templates are *not* normalized here
            logAge_grid[j, k] = logAge[j]
            metal_grid[j, k] = metal[k]

    return templates, lamRange_temp, logAge_grid, metal_grid, logLam2

#------------------------------------------------------------------------------

def ppxf_population(lams, data, error, z=0, inst='MUSE', sl_name='MILES_ppxf_default'):
    '''
    Master function that prepares the input to feed ppxf.
    
    :param: lams: input wavelength array
    :param: data: input spectrum in erg/s/cm**2/s/A
    :param: error: input error in (erg/s/cm**2/s/A)**2
    :param: inst: default='MUSE', what data is this ?
    :param: sl_name: which stellar libraries to use ?
    
    :return: ... lots of stuff ...
    '''
    
    lamRange = np.array([lams[0],lams[-1]])
    
    # Ok, I need to redshift this all back to the reference frame, or ppxf might crash.
    # Here's the original text from M. Cappellari:
    #
    # << If the galaxy is at a significant redshift (z > 0.03), one would need to apply
    # a large velocity shift in PPXF to match the template to the galaxy spectrum.
    # This would require a large initial value for the velocity (V > 1e4 km/s)
    # in the input parameter START = [V,sig]. This can cause PPXF to stop!
    # The solution consists of bringing the galaxy spectrum roughly to the
    # rest-frame wavelength, before calling PPXF. In practice there is no
    # need to modify the spectrum before the usual LOG_REBIN, given that a
    # red shift corresponds to a linear shift of the log-rebinned spectrum.
    # One just needs to compute the wavelength range in the rest-frame
    # and adjust the instrumental resolution of the galaxy observations. >>
    
    lams0 = lams/(z+1)
    lamRange0 = np.array([lams0[0],lams0[-1]])
    
    # Ok, now we need to make sure that we only fit the part of the spectrum that 
    # overlaps with the stellar libraries. 
    # We need to find the limits of these. Get one of them and figure it out.
    sl_fns = glob.glob(os.path.join(sl_models[sl_name]['sl_loc'],'*'))
    hdu_sl = pyfits.open(sl_fns[0])
    header_sl = hdu_sl[0].header
    hdu_sl.close()
    
    lams_sl = np.arange(0, header_sl['NAXIS1'],1)*header_sl['CDELT1'] + header_sl['CRVAL1']
    lamRange_sl = [lams_sl[0],lams_sl[-1]]
    
    fit_lims =[np.max([lamRange0[0],lamRange_sl[0]]), 
               np.min([lamRange0[1],lamRange_sl[1]])]
    
    # Make sure we only use the spectral region that is common to the data and the 
    # stellar libraries
    mask = (lams0 > fit_lims[0]) & (lams0 < fit_lims[1])
    data0 = data[mask]
    lams0 = lams0[mask]
    lamRange0 = np.array([lams0[0],lams0[-1]])
    
    # Need to log-rebin this ...
    galaxy, logLam0, velscale = util.log_rebin(lamRange0, data0)
    #galaxy = galaxy/np.median(galaxy)
    # I should probably do the same for the error eventually ...
    
    # The noise level is chosen to give Chi^2/DOF=1 without regularization (REGUL=0)
    #
    noise = galaxy*0 + 0.01528           # Assume constant noise per pixel here

	# Very well, now, let's get those templates ready
    templates, lamRange_temp, logAge_grid, metal_grid, logLam_temp = \
        setup_spectral_library(velscale, inst, sl_name, fit_lims, z)
	
	# From M. Cappellari:
    # << The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV.>>
    #
    dv = c*np.log(lamRange_temp[0]/lamRange0[0])  # km/s
    goodpixels = util.determine_goodpixels(logLam0, lamRange_temp, 0.0)
    
	# From M. Cappellari:
    # << IMPORTANT: Ideally one would like not to use any polynomial in the fit
    # as the continuum shape contains important information on the population.
    # Unfortunately this is often not feasible, due to small calibration
    # uncertainties in the spectral shape. To avoid affecting the line strength of
    # the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
    # multiplicative ones (MDEGREE=10). This is only recommended for population, not
    # for kinematic extraction, where additive polynomials are always recommended. >>
    #
    # Since I do not care too much about stellar kinematics (because my SN is usually low)
    # but worry a lot about my emission lines, I follow that advice here.
    
    vel = 0.   # Initial estimate of the galaxy velocity in km/s -it's been de-redshifted.
    start = [vel, 100.]  # (km/s), starting guess for [V,sigma]

    # See the pPXF documentation for the keyword REGUL,
    # for an explanation of the following two lines
    #
    templates /= np.median(templates) # Normalizes templates by a scalar
    regul_err = 0.004  # Desired regularization error

    t = clock()

    plt.clf()
    plt.subplot(211)

    pp = ppxf(templates, galaxy, noise, velscale, start,
              goodpixels=goodpixels, plot=False, moments=4, degree=-1,
              vsyst=dv, clean=False, mdegree=10, regul=1./regul_err)


    # When the two numbers below are the same, the solution is the smoothest
    # consistent with the observed spectrum.
    #
    '''
    print('Desired Delta Chi^2: %.4g' % np.sqrt(2*goodpixels.size))
    print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1)*goodpixels.size))
    print('Elapsed time in PPXF: %.2f s' % (clock() - t))

    print('Mass-weighted <logAge> [Gyr]: %.3g' %
          (np.sum(pp.weights*logAge_grid.ravel())/np.sum(pp.weights)))
    print('Mass-weighted <[M/H]>: %.3g' %
          (np.sum(pp.weights*metal_grid.ravel())/np.sum(pp.weights)))

    plt.subplot(212)
    s = templates.shape
    weights = pp.weights.reshape(s[1:])/pp.weights.sum()
    plt.imshow(np.rot90(weights), interpolation='nearest', 
               cmap='gist_heat', aspect='auto', origin='upper',
               extent=[np.log10(1), np.log10(17.7828), -1.9, 0.45])
    plt.colorbar()
    plt.title("Mass Fraction")
    plt.xlabel("log$_{10}$ Age (Gyr)")
    plt.ylabel("[M/H]")
    plt.tight_layout()
    plt.show()
    '''
    
    return logLam0,galaxy,pp
    
    
#------------------------------------------------------------------------------
