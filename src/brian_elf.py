# -*- coding: utf-8 -*-
#
# This file contains several function and tools used by the BRIAN routines to fit
# emission line.
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

from __future__ import print_function

import numpy as np
import sys
import os
import scipy.stats as stats
import cap_mpfit as mpfit

from brian_metadata import *


# ----------------------------------------------------------------------------------------

def inst_resolution(inst = 'MUSE', get_ff = False, show_plot=False):
    '''
    This function returns the resolution as a function of the wavelength
    for different instruments.
    
    Returns a callable function of the wavelength (in Angstroem !).
    '''
     
    if inst == 'MUSE':
       if get_ff:
          this_fn_path = os.path.dirname(__file__)
          ref_fn = 'MUSE_specs/MUSE_spectral_resolution.txt'
          R_lam = np.loadtxt(os.path.join(this_fn_path,ref_fn),
                           skiprows = 1)
                           
          # Fit a polynomial to this. Deg 3 works well.
          z = np.polyfit(R_lam[:,0]*10.,R_lam[:,1],3)
          p = np.poly1d(z)
          
          if show_plot:
             plt.close(99)
             plt.figure(99)
             
             lams = np.arange(4500,10000.,1)
             
             plt.plot(lams,p(lams), 'b-')            
             plt.plot(R_lam[:,0]*10.,R_lam[:,1], 'k.')
             
             plt.show()
          
       else:
           #Fit on 02.2016:
           p = np.poly1d([ -8.27037043e-09, 1.40175196e-04, -2.83940026e-01, 7.13549344e+02])
           
       return p    
                    
    else:
        sys.exit('Unknown instrument...')
# ----------------------------------------------------------------------------------------      

def obs_sigma(sigma, lam, inst='MUSE', z=0, in_errs = None):
    '''
    This function compute the observed sigma (in [A]) for an emission line, given a "real" 
    sigma in km/s and the line position. Set 'inst'=None to NOT account for the instrument
    dispersion.
    '''
    
    if in_errs:
        errs = in_errs
    else:
        errs = np.zeros(3)
    
    # Here, also account for the redshift. I.e. redshifted lines become a bit fatter.
    obj_sigma = sigma/c * lam * (1+z)
    # WARNING: THIS SHOULD ALSO BE CORRECTED FOR THE LINE REDSHIFT !
    # TODO IN THE REST OF CODE
    obj_sigma_err = obj_sigma**2 * (errs[0]/sigma**2 + 
                                        errs[1]/lam**2 + 
                                        errs[2]/(z+1)**2 )

    if inst:
        # TODO: understand what is the TRUE spectral resolution of MUSE !
        # Measure it from the sky emission lines + arc lines ?      
            
        # What is the velocity dispersion we have ? Convert from km/s to Angstroem
        # and account for the instrumental dispersion too !        
        inst_sigma = lam/(inst_resolution(inst=inst)(lam)*2*np.sqrt(2*np.log(2)))
        # Here, assume no error coming from the Resolution curve ... sigh...
        inst_sigma_err = errs[1]/(inst_resolution(inst=inst)(lam)*2*np.sqrt(2*np.log(2)))**2
        
        # Then combine sigma, in A
        this_sigma = np.sqrt(obj_sigma**2 + inst_sigma**2)
        this_sigma_err = inst_sigma_err * (inst_sigma/this_sigma)**2 + \
                         obj_sigma_err * (obj_sigma/this_sigma)**2
        
    else:
        this_sigma = obj_sigma
        this_sigma_err = obj_sigma_err


    if in_errs:        
        return [this_sigma,this_sigma_err]
    else:
        return this_sigma
# ----------------------------------------------------------------------------------------

def gauss_profile(x,I,mu,sigma):
    '''
    This is a normal gaussian, normalized to peak. mu and sigma in Angstroem
    '''
    
    y = I * stats.norm.pdf(x,loc=mu,scale=sigma)/stats.norm.pdf(0,loc=0,scale=sigma)
    # As usual, I normalize the function to peak, so that "I" corresponds to the line peak
    # intensity.
    
    return y 
# ----------------------------------------------------------------------------------------

def gauss_hist(x,binedges,I,mu,sigma):
    '''
    Compute the proper histogram of an emission line, given certain bin sizes. 
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
def els_spec(x, p, method ='gauss_hist', be=None, inst='MUSE'):
    '''
    Computes a pure emission line spectrum, given a set of parameters.
    'method' allows to choose between gaussian profile or gauss-hermite profile (not yet
    supported). 
    if gaussian profile: p = lam0 (no shift), I, v (km/s), sigma (km/s), etc ...
    if gauss-hermite :   p = lam0 (no shift), I, v (km/s), sigma (km/s), h3, h4 etc ...
    '''
    
    spec = np.zeros_like(x)
    if method=='gauss_hist':
        for i in range(len(p)/4):
            # What is the wavelength we are at ? Convert from km/s to Angstroem
            this_lam = (p[4*i+2]/c+1)*p[4*i]
            # Get the observed sigma in Angstroem, accounting for the instrumental res.
            this_sigma = obs_sigma(p[4*i+3],this_lam, inst=inst)
            spec += gauss_hist(x,be,p[4*i+1],this_lam,this_sigma)
            
    elif method=='gauss_profile':
        for i in range(len(p)/4):
            # What is the wavelength we are at ? Convert from km/s to Angstroem
            this_lam = (p[4*i+2]/c+1)*p[4*i]
            # Get the observed sigma in Angstroem, accounting for the instrumental res.
            this_sigma = obs_sigma(p[4*i+3],this_lam, inst=inst)
            spec += gauss_profile(x,be,p[4*i+1],this_lam,this_sigma)       
            
    #elif method =='gaussherm_hist':
    else:
        sys.exit('Method %s not yet supported !' % method)    
    
    return spec
# ----------------------------------------------------------------------------------------

def line_fit_erf(p,fjac=None, x=None,y=None,err=None, method='gauss_hist', be=None, 
                 inst='MUSE'):
    '''
    Computes the residuals to be minimized by mpfit, given a model and data.
    '''
    
    model = np.zeros_like(y)
    model = els_spec(x, p, method=method, be=be, inst=inst)    
    
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
    
    # Create the parinfo: each line has a ref wavelength, a peak intensity, a mean and 
    # a width. Note that the ref wavelength WILL NOT BE FITTED, IT IS FOR REF. ONLY !
    parinfo = [{'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} 
               for i in range(4*len(elines.keys()))]
    # And the list with our initial guesses
    p0 = [0 for i in range(len(elines.keys()*4))]

    # Calculate the references starting velocity
    ref_dv_line = params['ref_dv_line']

    # Extract a spectrum around that line
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
        parinfo[4*i]['fixed'] = 1
        p0[4*i] = elines[line][0][0]

        # Fill sigma and v inside p0
        p0[i*4+3] = elines[line][0][2]
        p0[i*4+2] = v0
	
        # For the intensity, pick the max of the area
        loc_spec = spec[(lams > elines[line][0][0] * (params['z_target'] + 1 -300./c))*
                        (lams < elines[line][0][0] * (params['z_target'] + 1 +300./c))]
        p0[i*4+1] = np.max(loc_spec)
    
        # And now fill all the associated constraints
        for (j,constraints) in enumerate(elines[line][1:]):
            for key in constraints.keys():
                # Decode the "tied" string, and re-format as required.
                if key =='tied' and constraints[key]:
                    ties = constraints[key].split('p(')
                    tied_line = np.where(np.sort(elines.keys()) == ties[-1][:-1])[0][0]
                    ties = ties[0]+'p['+np.str(4*tied_line+j+1)+']'
                    parinfo[4*i+j+1][key] = ties
                else:
                    # Assign the constraints
                    parinfo[4*i+j+1][key] = constraints[key]
        
    # Create the dictionary containing the spectra, error, etc ...    
    fa = {'x':lams, 'y':spec, 'err':np.sqrt(err), 'method':params['line_profile'], 
          'be':be, 'inst':params['inst']}

    m = mpfit.mpfit(line_fit_erf, p0, functkw=fa, parinfo=parinfo, quiet=1) 
    
    return m
# ----------------------------------------------------------------------------------------    