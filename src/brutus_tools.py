# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains general tools for the brutus routines to fit the stellar continuum and 
 the emission lines in an IFU data cube.

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import signal 

from brutus_metadata import c
from brutus_metadata import __version__ as version
 
 # ----------------------------------------------------------------------------------------      
  
def init_worker():
    '''Handles KeyboardInterrupt during multiprocessing.
    
    :Notes:
        See https://noswap.com/blog/python-multiprocessing-keyboardinterrupt
    '''
    signal.signal(signal.SIGINT, signal.SIG_IGN)

# ----------------------------------------------------------------------------------------      
  
def hdu_add_brutus(hdu,procstep):
    '''Adds dedicated brutus keywords to a FITS file header.
    
    :Args:
        hdu: FITS hdu
                The destination hdu to which the brutus keywords must be added.
        procstep: string
                  The name of the processing step creating the FITS file.
    
    :Returns:
        out: FITS header
             The newheader with brutus info included.   
    '''
    
    hdu.header['BRUTUS'] = (version, 'brutus version that created this file.')
    hdu.header['BRSTEP'] = (procstep, 'brutus processing step that created this file.')
    
    return hdu
    
# ----------------------------------------------------------------------------------------      
  
def hdu_add_wcs(newhdu,refheader):
    '''Adds the WCS coordinates from a reference header to a new hdu.
    
    :Args:
        newheader: FITS hdu
                   The destination hdu to which the WCS keywords must be added.
        refheader: FITS header
                   The reference header, from which to transer the WCS keywords.
    
    :Returns:
        out: FITS hdu
             The new hdu with WCS info included.   
    :Notes:
        Keywords transfered are  'CRPIX1', 'CD1_1', 'CTYPE1', 'CUNIT1', 'CRPIX2', 'CD2_2',
        'CTYPE2', 'CUNIT2', 'CD1_2', 'CD2_1', 'CRVAL1' and 'CRVAL2'.
    '''
    
    newhdu.header['CRPIX1'] = refheader['CRPIX1']
    newhdu.header['CD1_1'] = refheader['CD1_1']
    newhdu.header['CTYPE1'] = refheader['CTYPE1']
    newhdu.header['CUNIT1'] = refheader['CUNIT1']
    newhdu.header['CRPIX2'] = refheader['CRPIX2']
    newhdu.header['CD2_2'] = refheader['CD2_2']
    newhdu.header['CTYPE2'] = refheader['CTYPE2']
    newhdu.header['CUNIT2'] = refheader['CUNIT2']
    newhdu.header['CD1_2'] = refheader['CD1_2']
    newhdu.header['CD2_1'] = refheader['CD2_1']
    newhdu.header['CRVAL1'] = refheader['CRVAL1']
    newhdu.header['CRVAL2'] = refheader['CRVAL2']
    
    return newhdu
# ----------------------------------------------------------------------------------------      
    
def hdu_add_lams(newhdu,refheader):
    '''Adds the wavelength information from a reference header to a new hdu.
    
    :Args:
        newhdu: FITS hdu
                The destination hdu to which the wavelength keywords must be added.
        refheader: FITS header
                   The reference header, from which to transer the wavelength keywords.
    
    :Returns:
        out: FITS hdu
             The newheader with wavelength info included.   
    :Notes:
        Keywords transfered are 'CTYPE3', 'CUNIT3', 'CD3_3', 'CRPIX3', 'CRVAL3', 'CD1_3', 
        'CD2_3', 'CD3_1' and 'CD3_2'.
    '''
    
    newhdu.header['CTYPE3'] = refheader['CTYPE3'] 
    newhdu.header['CUNIT3'] = refheader['CUNIT3']   
    newhdu.header['CD3_3'] = refheader['CD3_3']
    newhdu.header['CRPIX3'] = refheader['CRPIX3']
    newhdu.header['CRVAL3'] = refheader['CRVAL3']
    newhdu.header['CD1_3'] = refheader['CD1_3']
    newhdu.header['CD2_3'] = refheader['CD2_3']
    newhdu.header['CD3_1'] = refheader['CD3_1']
    newhdu.header['CD3_2'] = refheader['CD3_2']
 
    return newhdu    
 
 # ----------------------------------------------------------------------------------------      
 
def spec_rebin(borders,new_borders,spec):
    ''' Rebin an array to new bins, while conserving the flux density.
    
    This function rebins a spectrum from one grid to another, assuming the two grids
    start and end at the same spot. Conserves the flux density in each bin.
    
    :Args:
        borders: 1-D numpy array
                 The edges of the initial spectrum bins.
        new_borders: 1-D numpy array
                     The edges of the new spectrum bins.
        spec: 1-D numpy array
              The spectrum to rebin, with len(spec) = len(borders)-1.
    
    :Returns:
        out: 1-D numpy array
             The rebinned spectrum, with size = len(new_borders)-1.
    '''
    if not( (borders[0]== new_borders[0]) * (borders[-1]== new_borders[-1])):
        sys.exit('Borders have different edges. Cannot rebin this.')
    
    # Find the indices where the new_border edges land in the old border edges. 
    # Not that with side="right", this gives us the "upper-bound of the bin the new_border
    # points land in.
    i_place = np.searchsorted(borders,new_borders, side='right').clip(0,len(borders)-1)
    
    # Create a storage structure for the new spectrum
    spec_new = np.zeros(len(new_borders)-1)
    
    # Start looping through the spectrum. There's a faster way (c.f. ppxf_util) but it
    # makes no sense to me. So here I take the long way, that I understand.
    for s in range (len(spec_new)):
        spec_new[s] = np.sum(spec[i_place[s]-1:i_place[s+1]] * \
                                            np.diff(borders)[i_place[s]-1:i_place[s+1]])
     
    # Correct for the first and last bin, which are not completely covered by the new bin
    spec_new -= (new_borders[:-1] - borders[i_place[:-1]-1]) * spec[i_place[:-1]-1]
    spec_new -= (borders[i_place[1:]] - new_borders[1:]) * spec[i_place[1:]-1]

    # Finally, normalize each bin. To maintain a uniform spectrum density
    spec_new /= np.diff(new_borders)
    
    return spec_new   
    
# ----------------------------------------------------------------------------------------      

def log_rebin (lams, spec, sampling = 1, velscale=None):
    '''Logarithmically rebin a spectrum.
    
    This function rebins a linear spectrum in log bins. It basically redistribute light
    into the new bins, while conserving the flux density. 
    
    :Args:
        lams: 1-D numpy array
              The input (linear) bin wavelengths, in Angstroem.
        spec: 1-D numpy array
              The input spectrum.
        sampling: int [default:1]
                  Amount of resampling (1 = no resampling).
        velscale: float [default: None]
                  The logarithmic bin size in km/s (None = choose default given sampling).
    
    :Returns:
        new_lams, new_spec, velscale: 1-D array, 1-D array, float
                                      The new spectral bin centers in Angstroem, 
                                      the reshaped spectrum, and the associated velscale 
                                      (in km/s).
    
    :Notes:
        Inspired by log_rebin inside ppxf_util, but without the one magic line I don't 
        understand (replaced by a "dumb" for-loop instead). Wraps around spec_rebin.  
        If velscale is set, sampling has no effect. The logarithmic bin centers are the 
        geometric means of the logarithmic bin edges. Unlike ppxf_util.log_rebin, return 
        the new bin centers in Angstroem.
    '''
    
    lam_range = [lams[0],lams[-1]]
    s = np.size(spec)
    
    dlam = (lam_range[1]-lam_range[0])/(s-1.) # Assumes constant lambda steps ...
    lim = lam_range/dlam + [-0.5,0.5] # All in units of dlam
    borders = np.linspace(lim[0],lim[1],s+1) # Current borders, in dlam steps
    
    loglim = np.log(lim)
    
    # No velscale selected. Return the default one.
    if velscale is None: 
        velscale = np.diff(loglim)/(sampling*s)*c   
        m = sampling*s
    else: # User has requested a certain velscale. Deal with it.
        logscale = velscale/c
        m = int(np.diff(loglim)/logscale)    # Number of output pixels
        loglim[1] = loglim[0] + m*logscale
    
    # Logarithmic borders
    new_borders = np.exp(np.linspace(*loglim, num=m+1)) 

    # Make sure the edges are the same, down to machine precision...
    new_borders[0] = borders[0]
    new_borders[-1] = borders[-1]
    
    # get the new spectrum, given those two grids
    new_spec = spec_rebin(borders,new_borders,spec)
    
    # Follow log_rebin in ppxf, and take the geometric average for the bin center
    new_lams = np.sqrt(new_borders[1:]*new_borders[:-1])*dlam
    
    return new_lams, new_spec, velscale
    
# ----------------------------------------------------------------------------------------      
    
def delog_rebin (lams, spec, old_lams, sampling = 1):
    '''Undoes the log_rebin function. 
    
    This function is designed to undo the brutus_tools.log_rebin function. Of course, this
    is not perfect, but setting sampling >1 helps avoid losses. 
    
    :Args:
        lams: 1-D numpy array
              The log-rebined wavelength array.
        spec: 1-D numpy array
              The log-rebinned spectra.
        old_lams: 1-D numpy array
              The original (linear) wavelength array.
        sampling: int [default: 1]
              Amount of resampling (1 = no resampling).
        
    :Returns:
        old_lams, old_spec: 1-D numpy array, 1-D numpy array
                            The delog-rebinned wavelength and spectra.

    :Notes:
        This function only deals with calculating the bin edges, and then feeds this to 
        spec_rebin.
    :See also:
        log_rebin()
    '''
    lam_range = [old_lams[0],old_lams[-1]]
    s = np.size(old_lams)
    
    dlam = (lam_range[1]-lam_range[0])/(s-1.) # Assumes constant lambda steps ...
    lim = lam_range/dlam + [-0.5,0.5] # All in units of dlam
    old_borders = np.linspace(lim[0],lim[1],s+1) # Current borders, in dlam steps
    
    loglim = np.log(lim)
    borders = np.exp(np.linspace(*loglim, num=sampling*s+1)) # Logarithmic borders

    # Make sure the edges are the same, down to machine precision...
    borders[0] = old_borders[0]
    borders[-1] = old_borders[-1]
    
    # I can simply get the new spectrum using the spec_rebin tool
    old_spec = spec_rebin(borders, old_borders, spec)
    
    return old_lams, old_spec
    
# ----------------------------------------------------------------------------------------

def inst_resolution(inst = 'MUSE', get_ff = False, show_plot=False):
    '''Returns the functional resolution of an instrument as a function of the wavelength.
    
    Returns a callable function of the wavelength (in Angstroem !).
    
    :Args:
        inst: string [default: 'MUSE']
              The name tag referring to a given instrument.
        get_ff: bool [default: False]
                Whether to recompute the given function from a reference dataset or not.
                Only valid with inst = 'MUSE'.
        show_plot: bool [default: False]
                   Whether to make a plot of the function.
    
    :Returns:
        R(lambda): function
                   A function that takes a float (lambda in Angstroem, and returns the 
                   corresponding value of the chosen instrument resolution.
    
    :Notes:
        Supported instruments: 'MUSE'
    
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
 
 
        