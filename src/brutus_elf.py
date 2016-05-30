# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains several function and tools used by the brutus routines to fit
 emission line.

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

from __future__ import print_function

import numpy as np
import scipy as sp
import warnings
import math
import sys
import os
import scipy.stats as stats
import brutus_mpfit as mpfit

from brutus_metadata import *
import brutus_tools

# ----------------------------------------------------------------------------------------      

def obs_sigma(sigma, lam, inst='MUSE', in_errs = None):
    '''
    This function compute the observed sigma (in [A]) for an emission line, given a "real" 
    sigma in km/s and the line position. Set 'inst'=None to NOT account for the instrument
    dispersion.
    '''
    
    # Do I want to propagate errors ?
    if in_errs:
        errs = in_errs
    else:
        errs = np.zeros(2)
     
    # From km/s to A    
    obj_sigma = sigma/c * lam
    # Note: here, lam is the redshifted wavelength. It includes the effect of
    # the redshift in fattening the line, i.e. sigma is the "true" sigma in the object
    # rest frame.
        
    obj_sigma_err = obj_sigma * np.sqrt(errs[0]**2/sigma**2 + errs[1]**2/lam**2)    
        
    if inst=='MUSE':
        
        # TODO: understand what is the TRUE spectral resolution of MUSE !
        # Measure it from the sky emission lines + arc lines ?      
        # For now, get it from the ref. file.    
        '''
        # What is the velocity dispersion we have in A?     
        inst_sigma = lam/(brutus_tools.inst_resolution(inst=inst)(lam) * 
                                                                   2*np.sqrt(2*np.log(2)))
        # Here, assume no error coming from the Resolution curve ... sigh...
        
        inst_sigma_err = errs[1]/(brutus_tools.inst_resolution(inst=inst)(lam) * 
                                                               2*np.sqrt(2*np.log(2)))**2
        
        # Then combine sigma, in A
        this_sigma = np.sqrt(obj_sigma**2 + inst_sigma**2)
        
        this_sigma_err = inst_sigma_err * (inst_sigma/this_sigma)**2 + \
                         obj_sigma_err * (obj_sigma/this_sigma)**2
        '''
        
        # For now, pretend that there is no instrumental dispersion
        this_sigma = obj_sigma
        this_sigma_err = obj_sigma_err
        
    else:
    
        # Pretend there is no instrumental dispersion
        this_sigma = obj_sigma
        this_sigma_err = obj_sigma_err

    if in_errs:        
        return [this_sigma,this_sigma_err]
    else:
        return this_sigma
# ----------------------------------------------------------------------------------------

def gauss_profile(x,I,mu,sigma):
    '''Computes a gaussian profile at x, normalized to peak.
    
    :Args:
        x: float, 1-D numpy array
           Locations to compute the Gaussian profile, assumed to be in Angstroem.
        I: float
           The profile peak intensity.
        mu: float
            The profile mean, in Angstroem.
        sigma: float
               The standard deviation of the profile, in Angstroem.
           
    :Returns:
        out: float, 1-D numpy array
             The profile evaluated at x.
    '''
    
    y = I * stats.norm.pdf(x,loc=mu,scale=sigma)/stats.norm.pdf(0,loc=0,scale=sigma)
    # As usual, I normalize the function to peak, so that "I" corresponds to the line peak
    # intensity.
    
    return y 
# ----------------------------------------------------------------------------------------

def gauss_hist(x,binedges,I,mu,sigma):
    '''Computes the proper histogram of an emission line, given certain bin sizes. 
    
    The profile is normalized to peak.
    
    :Args:
        x: float, 1-D numpy array
           Locations to compute the Gaussian profile, assumed to be in Angstroem.
        binedges: 1-D numpy array
           Locations of the bin edges, with size of len(x)+1.
        I: float
           The profile peak intensity.
        mu: float
            The profile mean, in Angstroem.
        sigma: float
               The standard deviation of the profile, in Angstroem.
           
    :Returns:
        out: float, 1-D numpy array
             The profile evaluated at x
             
    :Notes:
        This function is taking into account the finite bin sizes when observing a line 
        with an intrinsic gaussian profile. I.e. each location shows the mean intensity 
        over the bin, rather than the intensity of the gaussian profile at the bin center. 
        This is the functional form to use when fitting an emission line in a spectra, 
        rather than a simple gaussian profile.
    '''
    # First, use the fast cumulative distribution function (CDF) of the normal 
    # distribution to integrate a gaussian (normalized to peak!) between a set of bin 
    # edges, and normalize with the bin width. 
    # This is EXACTLY what the data in a flux-density spectrum should be: 
    # the average flux-density over the bin. 
    
    # 1) Integrate over the bins
    y = I * (stats.norm.cdf(binedges[1:],loc=mu,scale=sigma) - 
                 stats.norm.cdf(binedges[:-1],loc=mu, scale=sigma))
    # 2) Scale the distribution, to make sure "I" is the line peak value. stats.norm is by
    # default normalized, but we want 1 for x=0.
    y /= stats.norm.pdf(0,loc=0,scale=sigma)        
    # 2) Normalize over the bins sizes
    y /= np.diff(binedges)
    
    return y
# ----------------------------------------------------------------------------------------

def hermite_poly_scale(deg):
    '''
    Computes the scaling factor to be applied to Hermite polynomials, when using offcenter
    & broad gaussian line profiles (i.e. with mu !=0 and and sigma != 1).
    See footnote 3 on p. 555 of van der Marel (1993).
    '''
    return 1./(np.sqrt(math.factorial(deg)*2.**deg))
# ----------------------------------------------------------------------------------------

def gauss_hermite_profile(x,I,mu,sigma,h3,h4):
    ''' 
    Construct a Gauss-Hermite profile, following van der Marel (1993) and Riffel (2010)
    Set h2 = h1 = 0, h0 = 1 => it's a gaussian plus h3 and h4
    mu and sigma in Angstroem
    '''
    
    # First, the gaussian - normalized to peak.
    y = I * stats.norm.pdf(x,loc=mu,scale=sigma)/stats.norm.pdf(0,loc=0,scale=sigma)
	
    # Construct the Hermite polynomials - and scale the default scipy ones as required.
    H3 = hermite_poly_scale(3)*sp.special.eval_hermite(3,(x-mu)/sigma)
    H4 = hermite_poly_scale(4)*sp.special.eval_hermite(4,(x-mu)/sigma)
	
    y *= 1 + h3*H3 + h4*H4
	
    return y
# ----------------------------------------------------------------------------------------
	
def gauss_hermite_hist(x,binedges,I,mu,sigma,h3,h4):
    '''
    Computes the histogram profile of a line, given a Gauss-Hermite profile, and certain
    binedges for each sampling point.
    '''
    
    y = np.zeros_like(x)
    
    # Here, I need to integrate by hand, because I could not find a handy CDF function ...
    # Sigh, this is slow ... also because integrate.quad does not work with arrays ...
    # Re-Sigh ...
    for k in range(len(y)):
        y[k] = sp.integrate.quad(gauss_hermite_profile, # The function
                                         binedges[k],binedges[k+1], # The bin edges
                                         args=(I,mu,sigma,h3,h4) # The arguments
                                         )[0]
    
    y /= np.diff(binedges) # Normalize for each bin
    return y
# ----------------------------------------------------------------------------------------

# A function that constructs the spectrum of emission lines
def els_spec(x, p, method ='gauss', be=None, inst='MUSE'):
    '''
    Computes a pure emission line spectrum, given a set of parameters.
    'method' allows to choose between gaussian profile or gauss-hermite profile (not yet
    supported). 
    if gaussian profile: p = lam0 (no shift), I, v (km/s), sigma (km/s), etc ...
    if gauss-hermite :   p = lam0 (no shift), I, v (km/s), sigma (km/s), h3, h4 etc ...
    '''
    
    spec = np.zeros_like(x)
    nplp = 6 # Always use h3 and h4, even if they are fixed.
    
    for i in range(len(p)/nplp):
        # What is the wavelength we are at ? Convert from km/s to Angstroem
        this_lam = (p[nplp*i+2]/c+1.)*p[nplp*i]
        # Get the observed sigma in Angstroem, accounting for the instrumental res.
        this_sigma = obs_sigma(p[nplp*i+3],this_lam, inst=inst)
    
        if method=='gauss':
            spec += gauss_hist(x,be,p[nplp*i+1],this_lam,this_sigma)
            
        elif method=='gauss_profile': 
            # This is a "pure" gaussian profile. It does NOT take into account the fact that
            # we are measuring lines with a given bin -size. I.e. this assumes that the middle 
            # of the bin is at the "mean intensity" of the bin - which is not true for a 
            # gaussian, especially with crude bins !
        
            # This is for test purposes only - DON'T USE THIS UNLESS YOU KNOW WHAT THIS MEANS !
            spec += gauss_profile(x,p[nplp*i+1],this_lam,this_sigma)       
               
        elif method =='gauss-herm':
            spec += gauss_hermite_hist(x,be,p[nplp*i+1],this_lam,this_sigma,p[nplp*i+4],
                                                                            p[nplp*i+5])
    
        elif method =='gauss-herm_profile':
            # This is a "pure" gaussian-hermite profile. It does NOT take into account the fact 
            # that we are measuring lines with a given bin-size. I.e. this assumes that the 
            # middle of the bin is at the "mean intensity" of the bin - which is not true for  
            # a gaussian, especially with crude bins !
        
            # This is for test purposes only - DON'T USE THIS UNLESS YOU KNOW WHAT THIS MEANS !
            spec += gauss_hermite_profile(x,p[nplp*i+1],this_lam,this_sigma,p[nplp*i+4],
                                                                            p[nplp*i+5])
           
        else:
            sys.exit('Method %s unknown !' % method)    
    
    return spec
# ----------------------------------------------------------------------------------------

def line_fit_erf(p,fjac=None, x=None,y=None,err=None, method='gauss', be=None, 
                 inst='MUSE'):
    '''
    Computes the residuals to be minimized by mpfit, given a model and data.
    '''
    
    model = np.zeros_like(y)
    model = els_spec(x, p, method=method, be=be, inst=inst)    
    
    # TODO:
    # Can I speed things up by looking only at the area with emission lines ?
    
    status = 0
    
    return ([status,(y-model)/err])
# ----------------------------------------------------------------------------------------
    
def els_mpfit(specerr,lams=None, be=None, params=None):
    '''
    The master function calling mpfit. Its construction implicitly assumes a single 
    wavelength array over the entire cube !  
    '''
    
    # Because of the multiprocessing unit only able to feed a single variable, 
    # I combined both the spectra and error inside one list. Let's get them out.
    spec = specerr[0]
    err = specerr[1]

    # First of all, create the parinfo and p0 structures for the fit.
    elines = params['elines']
    
    # Assume a fixed size for the 'eline' paramt structure - always with h3 and h4, even 
    # if NOT used ...
    if params['line_profile'] in ['gauss','gauss_profile']:
        nplp = 6 
    elif params['line_profile'] in ['gauss-herm','gauss-herm_profile']:
        nplp = 6
    
    # Create the parinfo: each line has a ref wavelength, a peak intensity, a mean and 
    # a width. Note that the ref wavelength WILL NOT BE FITTED, IT IS FOR REF. ONLY !
    parinfo = [{'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
               for i in range(nplp * len(elines.keys()))]
               
    # And the list with our initial guesses
    p0 = [0 for i in range(len(elines.keys() * nplp))]

    # Calculate the references starting velocity
    ref_dv_line = params['ref_dv_line']

    # Extract a spectrum around that line
    # WARNING: this could be a problem, if the velocities in the data are larger than this
    loc_spec = spec[(lams > (ref_dv_line * (params['z_target'] + 1 -300./c))) * \
                    (lams < (ref_dv_line * (params['z_target'] + 1 +300./c)))]
    loc_lams = lams[(lams > (ref_dv_line * (params['z_target'] + 1 -300./c))) * \
                    (lams < (ref_dv_line * (params['z_target'] + 1 +300./c)))]
                    
    # Use the spec maximum with +- dv range for the start wavelength. 
    # Sensitive to bad pixels !
    # TODO: fit a more robust way of the initial dv guess !
    v0 = (loc_lams[np.argmax(loc_spec)]/ref_dv_line-1.) * c

    # Now loop through each line, and fill the parinfo structure according to what the 
    # user wants to do. Also fill the starting guesses. I sort the lines according to 
    # their keys. This is to ensure I keep the same order throughout.
    for (i,line) in enumerate(np.sort(elines.keys())):
        # Make sure the "reference wavelength of the line is fixed, and just carried over 
        # for "fun". This is not something to fit, but I need access to it later on.
        parinfo[nplp*i]['fixed'] = 1
        p0[nplp*i] = elines[line][0][0]

        # Fill sigma and v inside p0
        p0[i*nplp+3] = elines[line][0][2]
        p0[i*nplp+2] = v0
        
        # Make sure the starting guess for the velocity is within the boundaries
        if params['elines'][line][2]['limited'][0] ==1: # There is a lower boundary ...
            p0[i*nplp+2] = np.max([params['elines'][line][2]['limits'][0], p0[i*nplp+2]])
        if params['elines'][line][2]['limited'][1] ==1: # There is an upper boundary
            p0[i*nplp+2] = np.min([params['elines'][line][2]['limits'][1], p0[i*nplp+2]])
	
        # For the intensity, pick the max of the area
        this_spec = spec[(lams > elines[line][0][0] * (params['z_target'] + 1 -300./c))*
                             (lams < elines[line][0][0] * (params['z_target'] + 1 +300./c))]
        
        # Define the initial guess as the largest peak in the signal - account for the 
        # sign, to allow for both emission AND absorption lines !
        # I get some warnings for all-nans slices ... damn ... For clarity in the prompt, 
        # let's catch them and ignore them just this once, if the user is ok with it.
        with warnings.catch_warnings():
            warnings.simplefilter(params['warnings'], category=RuntimeWarning)
            this_extrema = np.nanmax(np.abs(this_spec))
            
        if np.isnan(this_extrema):
            # Bad spectra with only nans
            this_extrema= np.nan
        else:
            p0[i*nplp+1] = this_extrema * \
                            np.sign(this_spec[np.abs(this_spec)==this_extrema][0])
        
        # Make sure the starting guess for the intensity is within the boundaries
        if params['elines'][line][1]['limited'][0] ==1: # There is a lower boundary ...
            p0[i*nplp+1] = np.max([params['elines'][line][1]['limits'][0], p0[i*nplp+1]])
        if params['elines'][line][1]['limited'][1] ==1: # There is an upper boundary
            p0[i*nplp+1] = np.min([params['elines'][line][1]['limits'][1], p0[i*nplp+1]])
        

        # If gauss-hermite, let the h3 and h4 = 0 as starting guesses.
        # Make sure to fix those parameters to 0 if I only want to do a normal gauss fit !
        # Otherwise, do as the user pleases ...
        if not('herm' in params['line_profile']):
             parinfo[nplp*i+4]['fixed'] = 1
             parinfo[nplp*i+5]['fixed'] = 1
             p0[i*nplp+4] = 0.
             p0[i*nplp+5] = 0.
    
        # And now fill all the associated constraints
        for (j,constraints) in enumerate(elines[line][1:]):
            for key in constraints.keys():
                # Decode the "tied" string, and re-format as required.
                if key =='tied' and constraints[key]:
                    ties = constraints[key].split('p(')
                    tied_line = np.where(np.sort(elines.keys()) == ties[-1][:-1])[0][0]
                    ties = ties[0]+'p['+np.str(nplp*tied_line+j+1)+']'
                    parinfo[nplp*i+j+1][key] = ties
                else:
                    # Assign the constraints
                    parinfo[nplp*i+j+1][key] = constraints[key]
  
    # Create the dictionary containing the spectra, error, etc ...    
    fa = {'x':lams, 'y':spec, 'err':np.sqrt(err), 'method':params['line_profile'], 
          'be':be, 'inst':params['inst']}
          
    m = mpfit.mpfit(line_fit_erf, p0, functkw=fa, parinfo=parinfo, quiet=1) 

    return m
# ----------------------------------------------------------------------------------------    