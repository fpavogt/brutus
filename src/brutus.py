# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains the master brutus routines to fit the stellar continuum and the 
 emission lines in an IFU data cube (i.e. MUSE). Most of these routines call sub-routines,
 after setting the scene/loading datasets/etc ...
 
 Any processing step MUST have a dediacted routine in this file call 'run_XXX', which can
 then refer to any existing/new brutus/Python module.

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
 
 -----------------------------------------------------------------------------------------
  
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import sys
import os
from astropy.io import fits as pyfits
from functools import partial
import pickle
import multiprocessing
import warnings
import glob
import warnings

# Import brutus-specific tools
import brutus_tools
import brutus_cof
import brutus_elf
import brutus_plots
import brutus_red
from brutus_metadata import *
from brutus_metadata import __version__

# Import other important tools, but don't fail, in case the user doesn't need them.
try:
    import brutus_ppxf
    import ppxf_util as util
except:
    warnings.warn("Faild to load brutus_ppxf. Is ppxf installed ?")

try:
    import fit_kinematic_pa as fkp
except:
    warnings.warn("Faild to load fit_kinematic_pa. Is it installed ?")

try:
    import pyqz
except:
    warnings.warn("Faild to load pyqz. Is it installed ?")

# ---------------------------------------------------------------------------------------- 
  
def run_snr_maps(fn_list, params, suffix = None, do_plot = False):
    '''
    This function computes the SNR maps for the continuum and Ha (or other line). 
    It also creates a map of spaxels with any signal at all.
    The resulting maps are saved to a fits file with full header and WCS coordinates.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make a plot of the SNR maps or not.
        
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        This function is a "first guess" of the SNR for latter use in the code. A more
        reliable measure of the SNR for the emission line should be computed after they
        have been fitted.
    '''
    
    if params['verbose']:
        print '-> Computing the SNR maps.' 
    
    cont_range = params['cont_range']
    line_range = [params['ref_dv_line']*(1.-200./c), params['ref_dv_line']*(1.+200./c)]
    
    # Import the fits file 
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array - REST frame !
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    lams /= params['z_target']+1
    
    # Continuum: median intensity level across range vs std across range
    # I get some warnings for all-nans slices ... damn ... For clarity in the prompt, 
    # let's catch them and ignore them just this once, if the user is ok with it.
    with warnings.catch_warnings():
        warnings.simplefilter(params['warnings'], category=RuntimeWarning)
        # The signal
        cont_s = np.nanmedian(data[(lams>=cont_range[0])*(lams<=cont_range[1]),:,:],axis=0) 
        # The noise
        cont_n = np.nanstd(data[(lams>=cont_range[0])*(lams<=cont_range[1]),:,:],axis=0) 
    
    cont_snr = cont_s/cont_n
    # Also make sure this is always > 0
    cont_snr[cont_snr <0] = 0
    
    # I also want to compute the SNR map for the strongest emission line in there.
    # This line is defined via the 'ref_dv_line' in params (i.e. it is the same line used 
    # for initial v guess for the emission line fitting).
    # I get some warnings for all-nans slices ... damn ... For clarity in the prompt, 
    # let's catch them and ignore them just this once, if the user is ok with it.
    with warnings.catch_warnings():
        warnings.simplefilter(params['warnings'], category=RuntimeWarning)
        line_s = np.nanmax(data[ (lams>=line_range[0]) * (lams<=line_range[1]), :,:], axis=0)
    # For the line, I measure the line "peak" above the noise. This is NOT ideal, but it
    # "only" needs to ignore" pixels with no signals to save time during the fitting. 
    line_snr = line_s/cont_n
    # Make sure this is always >0
    line_snr[line_snr<0] = 0
    
    # And create a map with just spaxels that have any data (i.e. have been observed).
    anything = np.ones_like(data[0,:,:])
    anything[np.all(np.isnan(data),axis=0)] = np.nan
        
    # Very well, now let's create a fits file to save this as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(cont_snr)
    hdu2 = pyfits.ImageHDU(line_snr)
    hdu3 = pyfits.ImageHDU(anything)
    # Make sure the WCS coordinates are included as well
    for hdu in [hdu1,hdu2,hdu3]:
        hdu = brutus_tools.hdu_add_wcs(hdu,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu = brutus_tools.hdu_add_brutus(hdu,suffix)
    # For reference, also include the line/region this maps are based on
    hdu1.header['BRRANG'] = (np.str(cont_range), 'spectral range used for continuum SNR')
    hdu2.header['BRRANG'] = (np.str(line_range), 'spectral range used for line SNR') 
        
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1,hdu2,hdu3])
    fn_out = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_snr.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # And add the filename to the dictionary of filenames
    fn_list['snr_cube'] = suffix+'_'+params['target']+'_snr.fits'
                
    if do_plot:
        # Alright, let's take out the big guns ...            
        brutus_plots.make_2Dplot(fn_out,1, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_cont_snr.pdf'),
                            contours = [3,5,10,20], vmin=0,vmax=30,
                            cbticks=[0,3,5,10,20,30], 
                            cblabel = r'Continuum SNR %.2f\AA\ - %.2f\AA' % (cont_range[0],cont_range[1]),
                            scalebar = params['scalebar'],
                            )  
        brutus_plots.make_2Dplot(fn_out,2, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_line_snr.pdf'),
                            contours = [3,5,10,20], vmin=0,vmax=30,
                            cbticks=[0,3,5,10,20,30],
                            cblabel = r'%.2f\AA\ emission line SNR' % params['ref_dv_line'],
                            scalebar = params['scalebar'],
                            )    
        brutus_plots.make_2Dplot(fn_out,3, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_signal.pdf'),
                            vmin=0,vmax=1,cbticks=[0,1],
                            cmap ='alligator',
                            cblabel = r'Spaxels with data',
                            scalebar = None
                            )               
                         
    return fn_list
# ---------------------------------------------------------------------------------------- 
  
def run_gal_dered(fn_list, params, suffix = None, do_plot = False):
    '''Corrects for Galactic extinction, given the Ab and Av extinction.
    
    This function erives the Alambda value for any wavelength, and corrects the data to
    correct for our "local" extinction. Intended for extragalactic sources.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make a plot of the Alambda correction applied or not.
        
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        To reproduce the approach from NED, use the Ab and Av value for you object from 
        there, and set curve='f99', rv=3.1.   
    '''
    
    if params['verbose']:
        print '-> Correcting for Galactic extinction.' 
    
    # Import the raw data 
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
    # Compute Alambda
    alams = brutus_red.alam(lams, params['Ab'],params['Av'], curve=params['gal_curve'],
                           rv=params['gal_rv'])
    
    # Compute the flux correction factor
    etau = brutus_red.galactic_red(lams,params['Ab'],params['Av'], 
                                  curve=params['gal_curve'],
                                  rv=params['gal_rv'])
    
    # Correct the data and the propagate the errors
    data *= etau[:, np.newaxis, np.newaxis]
    error *= (etau[:, np.newaxis, np.newaxis])**2
    
    # Save the datacubes
    hdu0 = pyfits.PrimaryHDU(None,header0)
    
    hdu1 = pyfits.ImageHDU(data)
    hdu2 = pyfits.ImageHDU(error)
        
    for hdu in [hdu1,hdu2]:
        # Make sure the WCS coordinates are included as well
        hdu = brutus_tools.hdu_add_wcs(hdu,header1)
        hdu = brutus_tools.hdu_add_lams(hdu,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu = brutus_tools.hdu_add_brutus(hdu,suffix)
        # Add the line reference wavelength for future references
        hdu.header['BR_AB'] = (params['Ab'], 'Extinction in B')
        hdu.header['BR_AV'] = (params['Av'], 'Extinction in V')
            
            
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1,hdu2])
    fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_galdered.fits')
    hdu.writeto(fn_out, clobber=True)
        
    # Add the filename to the dictionary of filenames
    fn_list['galdered_cube'] = suffix+'_'+params['target']+'_galdered.fits'
    
    # Do I want to save a plot ?
    if do_plot:
        ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                            '_gal_Alambda_corr.pdf')
        brutus_plots.make_galred_plot(lams,alams,etau,params['Ab'],params['Av'], ofn)
    
    return fn_list
    
# ----------------------------------------------------------------------------------------

def run_fit_continuum(fn_list, params, suffix=None, start_row = None, end_row = None, 
                      method='lowess'):
    ''' 
    This function fits the continuum in the datacube, either using ppxf (if SNR is decent)
    or using a simple polynomial value. It is designed to use multiprocessing to speed
    things up on good computers. It deals with the data columns-per-columns, and 
    can be restarted mid-course, in case of a crash. 
    '''
    
    # For the bad continuum: run a LOWESS filter from statsmodel.nonparametric
    #sm.nonparametric.smoothers_lowess.lowess(spec,lams,frac=0.05, it=5)
    # For the good fits, run ppxf
    
    # Rather than launch it all at once, let's be smart in case of problems. I'll run
	# the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
	# between each row.
	
	# First, load the datacube to be fitted. Use the galactic deredened one, if it exists.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
    # I also need to load the SNR cube for the spaxel selection
    hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['snr_cube']))
    snr_cont = hdu[1].data
    hdu.close()

    # Get some info about the cube
    nrows = header1['NAXIS1']
    if start_row is None:
        start_row = 0
    if end_row is None:
	    end_row = nrows-1
	
	# Ok, what do I want to do ?
    if method == 'lowess':
        if params['verbose']:
            print '-> Starting the continuum fitting using the LOWESS approach.'
            
        fit_func = partial(brutus_cof.lowess_fit, lams=lams, 
                           frac=params['lowess_frac'], it=params['lowess_it'])
        # Note here the clever use of the partial function, that turns the lowess_fit
        # function from something that takes 4 arguments into something that only takes 1
        # argument ... thus perfect for the upcoming "map" functions !
        
    elif method == 'ppxf':
        if params['verbose']:
            print '-> Starting the continuum fitting using PPXF.' 
        	
        # I need to do some special preparation for 'ppxf'. Namely, log-rebin the 
        # templates. I could do that for each fit, bit it would cost me time. At the price 
        # of a larger code here, do it only once, and save "tons" of time (~0.5s per fit).
        # The code below is based on the examples provided by M. Cappellari within ppxf 
        # istelf (ppxf_kinematics_example_sauron.py & ppxf_population_example_sdss.py), 
        # but has been modified to account a wavelength dependant spectral resolution of 
        # the data. And avoid doing the same thing multiple times via multiprocessing. 
        
        # The full spectral range
        lam_range = np.array([lams[0],lams[-1]])
        
        # De-redshift it using the user's best guess
        # Here's the original text from M. Cappellari:
        #
        # << If the galaxy is at a significant redshift (z > 0.03), one would need to 
        # apply a large velocity shift in PPXF to match the template to the galaxy 
        # spectrum. This would require a large initial value for the velocity (V>1e4 km/s)
        # in the input parameter START = [V,sig]. This can cause PPXF to stop!
        # The solution consists of bringing the galaxy spectrum roughly to the
        # rest-frame wavelength, before calling PPXF. In practice there is no
        # need to modify the spectrum before the usual LOG_REBIN, given that a
        # red shift corresponds to a linear shift of the log-rebinned spectrum.
        # One just needs to compute the wavelength range in the rest-frame
        # and adjust the instrumental resolution of the galaxy observations. >>
        lams0 = lams/(params['z_target']+1)
        lam_range0 = np.array([lams0[0],lams0[-1]])
        
        # We can only fit the spectra where they overlap with the spectral library.
        # Get one of the templates, figure out its range, and derive the fit limits.
        sl_fns = glob.glob(os.path.join(sl_models[params['ppxf_sl_name']]['sl_loc'],'*'))
        hdu_sl = pyfits.open(sl_fns[0])
        header_sl = hdu_sl[0].header
        hdu_sl.close()
        
        lams_sl = np.arange(0, header_sl['NAXIS1'],1)*header_sl['CDELT1'] + \
                                                                       header_sl['CRVAL1']
        lam_range_sl = [lams_sl[0],lams_sl[-1]]
        
        # What are my fit limits, then ?
        fit_lims =[np.max([lam_range0[0],lam_range_sl[0]]), 
                   np.min([lam_range0[1],lam_range_sl[1]])]  
        mask = (lams0 > fit_lims[0]) * (lams0 < fit_lims[1])
        lams0c = lams0[mask]
        lam_range0c = np.array([lams0c[0], lams0c[-1]])
        
        # Soon, I will be log-rebining the spectra. What are the new bins going to be ?
        log_lams0c, emptyspec, velscale = brutus_tools.log_rebin(lams0c,np.zeros_like(lams0c), 
                                                     sampling = params['ppxf_sampling'])
        
        # Deal with the template spectra: disperse them according to the instrument, and
        # log-rebin them
        if params['verbose']:
            sys.stdout.write('\r   Preparing the templates ... ')
            sys.stdout.flush()
            
        templates, lam_range_temp, logAge_grid, metal_grid, log_lams_temp = \
        brutus_ppxf.setup_spectral_library(velscale, params['inst'], 
                                          params['ppxf_sl_name'], fit_lims, 
                                          params['z_target'])
        if params['verbose']:   
            sys.stdout.write('done.')
            print ' '
            sys.stdout.flush()                                      
                                          
                                          
        # For the fit reconstruction later on, save the various wavelength ranges, etc ...
        fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_ppxf_lams.pkl')
        file = open(fn,'w')
        pickle.dump([lams,lams0,lams0c,log_lams0c],file)
        file.close()
        # And add the generic pickle filename to the dictionary of filenames
        fn_list['ppxf_lams'] = suffix+'_'+params['target']+'_ppxf_lams.pkl'              
                                          
        # From M. Cappellari:
        # << The galaxy and the template spectra do not have the same starting wavelength.
        # For this reason an extra velocity shift DV has to be applied to the template
        # to fit the galaxy spectrum. We remove this artificial shift by using the
        # keyword VSYST in the call to PPXF below, so that all velocities are
        # measured with respect to DV.>>
        #
        dv = c*np.log(lam_range_temp[0]/lam_range0c[0]) 
        
        # Initial estimate of the galaxy velocity in km/s:        
        vel = 0.   # It's been de-redshifted (using user's guess)
        start = [vel, 100.]  # (km/s), starting guess for [V,sigma]
        
        # Now create a list of "good pixels", i.e. mask emission lines and stuff.
        # List is the default one from ppxf.
        goodpixels = util.determine_goodpixels(np.log(log_lams0c), lam_range_temp, vel/c)

        # See the pPXF documentation for the keyword REGUL, and an explanation of the 
        # following two lines
        templates /= np.median(templates) # Normalizes templates by a scalar
        regul_err = params['ppxf_regul_err']  # Desired regularization error
        
        fit_func = partial(brutus_ppxf.ppxf_population, templates=templates, 
                           velscale=velscale, start=start, goodpixels=goodpixels, 
                           plot=False, moments=params['ppxf_moments'], 
                           degree=params['ppxf_degree'], vsyst=dv, clean=False, 
                           mdegree=params['ppxf_mdegree'], regul=1./regul_err)    

    # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
    # loose the current row. Just live with it.
    for row in np.linspace(start_row,end_row, end_row-start_row+1):   
	
		# Alright, now deal with the spaxels outside the user-chosen SNR range.
		# Replace them with nan's
        good_spaxels = np.ones((header1['NAXIS2']))
        if params[method+'_snr_min']:
            good_spaxels[snr_cont[:,row] < params[method+'_snr_min']] = np.nan
        if params[method+'_snr_max']:
            good_spaxels[snr_cont[:,row] >= params[method+'_snr_max']] = np.nan
		
		# Build a list of spectra to be fitted
        specs = [data[:,i,row] * good_spaxels[i] for i in range(header1['NAXIS2'])]
        errs = [error[:,i,row] * good_spaxels[i] for i in range(header1['NAXIS2'])]
        
        if method == 'ppxf':
            if params['verbose']:
                sys.stdout.write('\r   Log-rebin spectra in row %2.i ...            '%row)
                sys.stdout.flush()
                
            # I need to do some special preparation for 'ppxf'. Namely, log-rebin the 
            # spectra, and crop them to an appropriate wavelength after de-redshifting. 
            # Let's get started. 
            # To save time, only log-rebin spectra that are not all nans
            specs = [brutus_tools.log_rebin(lam_range0c, this_spec[mask],
                                     sampling=params['ppxf_sampling'])[1] 
                                     if not(np.all(np.isnan(this_spec))) else np.nan 
                                     for this_spec in specs]                      
            # Also take care of the error
            errs = [brutus_tools.log_rebin(lam_range0c, this_err[mask], 
                                     sampling=params['ppxf_sampling'])[1]
                                     if not(np.all(np.isnan(this_err))) else np.nan  
                                     for this_err in errs]                                                                 
                                                                            
            # Combine both spectra and errors, to feed only 1 element to fit_func
            # Turn the error into a std to feed ppxf
            specs = [[specs[i],np.sqrt(errs[i])] for i in range(header1['NAXIS2'])]
        
		# Set up the multiprocessing pool of workers
        if params['multiprocessing']:
            # Did the user specify a number of processes to use ?
            if type(params['multiprocessing']) == np.int:
                nproc = params['multiprocessing']
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
                
			# Limiting the maximum number of cpus for ppxf, to avoid memory errors... 
            if nproc > 8 and method == 'ppxf':
                nproc = 8
                
            pool = multiprocessing.Pool(processes = nproc, 
			                            initializer = brutus_tools.init_worker())    
			
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
                sys.stdout.flush()
			
            # Launch the fitting ! Make sure to deal with KeyBoard Interrupt properly
			# Only a problem for multiprocessing. For the rest of the code, whatever.
            try:   
                conts = pool.map(fit_func, specs)
            except KeyboardInterrupt:
                print ' interrupted !'
                # Still close and join properly
                pool.close()
                pool.join()
                sys.exit('Multiprocessing continuum fitting interrupted at row %i'% row)
            else: # If all is fine
                pool.close()
                pool.join()  
              
              
        else: # just do things 1-by-1
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, one at a time ...' % 
                                 row)
                sys.stdout.flush()
                                
            conts = map(fit_func, specs)
	    
	    # Here, I need to save these results. Pickle could be fast and temporary,
	    # Until I then re-build the entire cube later on ? Also allow for better
	    # row-by-row flexibility.
        fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_'+method+'_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')
        file = open(fn,'w')
        pickle.dump(conts,file)
        file.close()
        
    # And add the generic pickle filename to the dictionary of filenames
    fn_list[method+'_pickle'] = suffix+'_'+params['target']+'_'+method+'_row_'
   	     
    print ' done !'
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_make_continuum_cube(fn_list, params, suffix=None, method='lowess'):   
    ''' 
    This function is designed to construct a "usable and decent" datacube out of the
    mess generated by the continuum fitting function, i.e. out of the many pickle 
    files generated.
    '''
    
    if params['verbose']:
        print '-> Constructing the datacube for the continuum fitting (%s).' % method
       
    # First, load the original datacube. I need to know how much stuff was fitted.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    nrows = header1['NAXIS1']
    
    cont_cube = np.zeros_like(data) * np.nan

    # In the case of ppxf, also extract the stellar moments.
    # Not the most trustworthy, but interesting nonetheless.
    if params['ppxf_moments']>0:
        ppxf_sol_map = np.zeros((params['ppxf_moments'],header1['NAXIS2'],
                                                        header1['NAXIS1'])) * np.nan
    
    # For ppxf, also create a spectrum with the de-logged rebin original spectra. Useful 
    # to characterize the error associated with log-bin-delog-rebin process employed here.
    delog_raw_cube = np.zeros_like(data) * np.nan
                                                        
    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cube [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          fn_list[method+'_pickle']+str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            conts = pickle.load(myfile)  
            myfile.close() 
        
            # Mind the shape
            if method=='lowess':
                # I directly saved the continuum as an array just get it out.
                cont_cube[:,:,row] = np.array(conts).T 
            elif method =='ppxf':
                # Alright, in this case, the pickle file contains the output from ppxf
            
                # Extract the various wavelength ranges, etc ...
                fn = os.path.join(params['tmp_loc'],fn_list['ppxf_lams'])
                f = open(fn,'r')
                [lams,lams0,lams0c,log_lams0c] = pickle.load(f)
                f.close()
                
                # The fit outputs are all on the log_lams0c grid, and I want them back on
                # lams.  
                # First, extract the bestfit spectra. 
                # For test purposes, also extract the galaxy spectra with pp.galaxy
                bestfits = [pp[0] if not(pp is None) else None for pp in conts]
                galaxys = [pp[1] if not(pp is None) else None for pp in conts]
                
                # And the stellar moments
                ppxf_sol_map[:,:,row] = np.array([pp[2] if not(pp is None) else 
                                                  np.zeros(params['ppxf_moments'])*np.nan 
                                                  for pp in conts]).T
                                                  
                # Now, I need storage that spans the original data range
                full_bestfits = [np.zeros_like(lams0) for this_spec in bestfits]
                full_galaxys = [np.zeros_like(lams0) for this_spec in galaxys]
                
                # Now, I want to delog-rebin the spectra back to the original grid
                for (k,cube) in enumerate([bestfits, galaxys]):
                    cube = [brutus_tools.delog_rebin(log_lams0c, this_spec, lams0c, 
                                                    sampling=params['ppxf_sampling'])[1] 
                            if not(this_spec is None) else None for this_spec in cube]
                            
                    # Now, I want to un-crop this, by adding nans elsewhere. Make a loop, 
                    # but there's probably a more Pythonic way of doing it. 
                    # TODO: try again when my brain won't be that fried ...
                    mask = (lams0 >=lams0c[0]) * (lams0<=lams0c[-1])
                
                    for (f,fit) in enumerate(cube):
                        if not(fit is None):
                            [full_bestfits,full_galaxys][k][f][mask] = fit
                
                    # And finally, I need to de-redshift that spectra. No need to touch 
                    # the spectra because I did not touch it before.
                
                    # Ready to save the continuum cube
                    [cont_cube,delog_raw_cube][k][:,:,row] = \
                                         np.array([full_bestfits, full_galaxys][k]).T
    
    if method == 'ppxf':                                     
        # For consistency with the rest of the code, add the guess redshift to the velocity. 
        ppxf_sol_map[0,:,:] += params['z_target']*c
                                     
    # Very well, now let's create a fits file to save this as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(cont_cube)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
    hdu1 = brutus_tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_'+method+'.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # And add the filename to the dictionary of filenames
    fn_list[method+'_cube'] = suffix+'_'+params['target']+'_'+method+'.fits'
    
    if method =='ppxf':
    
        # Also add the delog-rebin raw cube, and the moments maps
        hdu0 = pyfits.PrimaryHDU(None,header0)
        hdu1 = pyfits.ImageHDU(delog_raw_cube)
        # Make sure the WCS coordinates are included as well
        hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
        hdu1 = brutus_tools.hdu_add_lams(hdu1,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
        hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_'+method+'_delog_raw.fits')
        hdu.writeto(fn_out, clobber=True)
    
        # And add the filename to the dictionary of filenames
        fn_list['ppxf_delog_raw_cube'] = suffix+'_'+params['target']+'_'+method+\
                                           '_delog_raw.fits'                     
        # ---                                   
        hdu0 = pyfits.PrimaryHDU(None,header0)
        hdu1 = pyfits.ImageHDU(ppxf_sol_map)
        # Make sure the WCS coordinates are included as well
        hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
        hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_'+method+'_sol_map.fits')
        hdu.writeto(fn_out, clobber=True)
    
        # And add the filename to the dictionary of filenames
        fn_list['ppxf_sol_map'] = suffix+'_'+params['target']+'_'+method+\
                                           '_sol_map.fits'                                  
                                           
    print ' '
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_plot_ppxf_sol(fn_list, params, suffix=None, vrange=None, sigrange=None):   
    ''' Creates some plots for the stellar kinematics parameters derived via ppxf.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        vrange: list of int [default: None]
                If set, the range of the colorbar for the velocity plot.
        sigrangre: list of int [default: None]        
                If set, the range of the colorbar for the velocity dispersion plot.
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames.
    '''    
    if params['verbose']:
        print '-> Making some nifty plots to visualize the output of ppxf.'
     
    # Open the file with the ppxf parameters   
    fn = os.path.join(params['prod_loc'], fn_list['ppxf_sol_map'])
    
    # Open the file
    hdu = pyfits.open(fn)
    header0 = hdu[0].header
    header1 = hdu[1].header
    data = hdu[1].data
    hdu.close()
    
    # Create a temporary FITS file
    fn_tmp = os.path.join(params['tmp_loc'],suffix+'_tmp.fits')

    # Because of the way aplpy works, I need to stored each "image" in its own fits
    # file. I don't want to keep them, so let's just use a temporary one, get the plot
    # done, and remove it. Not ideal, but who cares ?
        
    # Make single pretty plots for v, sigma, h3,h4, h5, h6 if they exist  
    type = ['v', 'sigma', 'h3', 'h4', 'h5', 'h6' ]                                     
    for t in range(np.shape(data)[0]):
        # Create a dedicated HDU
        tmphdu = pyfits.PrimaryHDU(data[t])
        # Add the WCS information
        tmphdu = brutus_tools.hdu_add_wcs(tmphdu,header1)
        tmphdu.writeto(fn_tmp, clobber=True)
            
        # Now plot it 
        this_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+'_ppxf_'+
                                                        type[t]+'.pdf') 
            
        # For the velocity fields, set the vmin and vmax
        if t == 0 and vrange:
            my_vmin = vrange[0]
            my_vmax = vrange[1]
            my_cmap = 'alligator'
            my_stretch = 'linear'
            my_label = r'$v$ [km s$^{-1}$]'
            my_cbticks = None
        elif t ==1 and sigrange:
            my_vmin = sigrange[0]
            my_vmax = sigrange[1]
            my_cmap = 'alligator'
            my_stretch = 'linear'
            my_label = r'$\sigma_{tot}$ [km s$^{-1}$]'
            my_cbticks = None
        else:
            my_vmin = None
            my_vmax = None
            my_cmap = None
            my_stretch = 'linear'
            my_label = ''
            my_cbticks = None
                                                            
        brutus_plots.make_2Dplot(fn_tmp,ext=0, ofn=this_ofn, contours=False, 
                                    vmin=my_vmin, vmax = my_vmax, cmap=my_cmap,
                                    stretch = my_stretch, cblabel=my_label, 
                                    cbticks = my_cbticks,
                                    scalebar = params['scalebar'],
                                    )                                   
    
    # Delete the temporary fits file
    os.remove(fn_tmp)
    
    return fn_list    
# ----------------------------------------------------------------------------------------

def run_build_continuum_mix(fn_list, params, suffix=None):   
    ''' Construct the mixed continuum cube, following the user's wishes.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames.
    '''   
    
    if params['verbose']:
        print '-> Constructing the mixed continuum cube.'
    
    # First, the raw datacube. Use the galactic deredened one if it exists.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # I also need the continuum cubes
    if fn_list['lowess_cube']:
        fn = os.path.join(params['prod_loc'],fn_list['lowess_cube'])
        if os.path.isfile(fn):
            hdu = pyfits.open(fn)
            cont_lowess = hdu[1].data
            hdu.close()
    
    if fn_list['ppxf_cube']:
        fn = os.path.join(params['prod_loc'],fn_list['ppxf_cube'])
        if os.path.isfile(fn):
            hdu = pyfits.open(fn)
            cont_ppxf = hdu[1].data
            hdu.close()
    
    # I also need to load the SNR cube for the spaxel selection
    hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['snr_cube']))
    snr_cont = hdu[1].data                               
    snr_elines = hdu[2].data
    hdu.close()        
    
     
    # Here, I don't use nan's. Should I ? 
    cont_mix = np.zeros_like(data) 
    
    # Very well, now, combine the different continua
    for key in params['which_cont_sub'].keys():
        if params['verbose']:
            print '   SNR %s: %s' % (key,params['which_cont_sub'][key])
            sys.stdout.flush()
                
        llim = np.int(key.split('->')[0])
        if key.split('->')[1] == 'max':
            ulim = np.nanmax(snr_cont)+1.
        else:
        	ulim = np.int(key.split('->')[1])
         
        # Assume the continuum subtraction is error free (either lowess or models).
        if params['which_cont_sub'][key] == 'lowess':
            cont_mix[:,(snr_cont>=llim)*(snr_cont<ulim)] = \
                     cont_lowess[:,(snr_cont>=llim) * (snr_cont<ulim)] 
        elif params['which_cont_sub'][key] == 'ppxf':
            cont_mix[:,(snr_cont>=llim)*(snr_cont<ulim)] = \
                     cont_ppxf[:,(snr_cont>=llim) * (snr_cont<ulim)]   
    
    # Save the cube
    hdu0 = pyfits.PrimaryHDU(None, header0)
    hdu1 = pyfits.ImageHDU(cont_mix)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
    hdu1 = brutus_tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_mixed_continuum_cube.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # Add the generic filename to the dictionary of filenames
    fn_list['cont_mix_cube'] = suffix+'_'+params['target']+'_mixed_continuum_cube.fits'
    
    return fn_list  
# ----------------------------------------------------------------------------------------

def run_fit_elines(fn_list, params, suffix=None, start_row = None, end_row = None, 
                   ):
    ''' 
    This function fits the emission lines in the datacube, after subtracting the continuum
    derived using LOWESSS or PPXF. It is designed to use multiprocessing to speed
    things up on good computers. It deals with the data columns-per-columns, and 
    can be restarted mid-course, in case of a crash. 
    '''
    
    nlines = len(params['elines'].keys())
    if params['verbose']:
        print '-> Starting the emission line fitting for %2.i line(s).' % nlines 
    
    # Rather than launch it all at once, let's be smart in case of problems. I'll run
	# the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
	# between each row.
	
	# First, load the datacube to be fitted. Use the galactic deredened one if it exists.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    # Calculate the outer bin edges for the spectrum
    be = np.append(lams-header1['CD3_3']/2.,lams[-1]+header1['CD3_3']/2.)
    
    # I also need to load the SNR cube for the spaxel selection
    hdu = pyfits.open(os.path.join(params['prod_loc'],
                                   fn_list['snr_cube']))
    snr_cont = hdu[1].data                               
    snr_elines = hdu[2].data
    hdu.close()
    
    # I also need the continuum cubes
    if fn_list['cont_mix_cube']:
        fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
        if os.path.isfile(fn):
            hdu = pyfits.open(fn)
            cont_mix_cube = hdu[1].data
            hdu.close()
        
        # Very well, now, perform the continuum subtraction.
        data -= cont_mix_cube
    
    # Get some info about the cube
    nrows = header1['NAXIS1']
    if start_row is None:
        start_row = 0
    if end_row is None:
	    end_row = nrows-1 
        
    fit_func = partial(brutus_elf.els_mpfit, lams=lams, be=be, params=params)
	# Note here the clever use of the partial function, that turns the els_mpfit
	# function from something that takes 4 arguments into something that only takes 1
	# argument ... thus perfect for the upcoming "map" functions !
    
    # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
	# loose the current row. Just live with it.
    for row in np.linspace(start_row,end_row, end_row-start_row+1):   
	    
		# Alright, now deal with the spaxels outside the user-chosen SNR range.
		# Replace them with nan's
        good_spaxels = np.ones((header1['NAXIS2']))
        if params['elines_snr_min']:
             good_spaxels[snr_elines[:,row]<params['elines_snr_min']] = np.nan
        if params['elines_snr_max']:
             good_spaxels[snr_elines[:,row]>params['elines_snr_max']] = np.nan
		
		# Build a list of spectra to be fitted
        specerrs = [[data[:,i,row] * good_spaxels[i],
                     error[:,i,row]* good_spaxels[i]] for i in range(header1['NAXIS2'])]
        
		# Set up the multiprocessing pool of workers
        if params['multiprocessing']:
            # Did the user specify a number of processes to use ?
            if type(params['multiprocessing']) == np.int:
                nproc = params['multiprocessing']
                pool = multiprocessing.Pool(processes = nproc, 
                                            initializer = brutus_tools.init_worker())
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(processes=None, 
                                            initializer = brutus_tools.init_worker())
			
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
                sys.stdout.flush()
			
			# Launch the fitting ! Make sure to deal with KeyBoard Interrupt properly
			# Only a problem for multiprocessing. For the rest of the code, whatever.
            try:   
                els = pool.map(fit_func, specerrs)
            except KeyboardInterrupt:
                print ' interrupted !'
                # Still close and join properly
                pool.close()
                pool.join()
                sys.exit('Multiprocessing line fitting interrupted at row %i'% row)
            else: # If all is fine
                pool.close()
                pool.join()  
              
        else: # just do things 1-by-1
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, one at a time ...' % 
                                 row)
                sys.stdout.flush() 
                                
            els = map(fit_func, specerrs)
	    
	    # Here, I need to save these results. Pickle could be fast and temporary,
	    # Until I then re-build the entire cube later on ? Also allow for better
	    # row-by-row flexibility.
        fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_elines_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')
        file = open(fn,'w')
        pickle.dump(els,file)
        file.close()
   	
   	# Add the generic filename to the dictionary of filenames
   	fn_list['elines_pickle'] = suffix+'_'+params['target']+'_elines_row_'
   	     
    print ' done !'
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_make_elines_cube(fn_list, params, suffix=None):   
    ''' 
    This function is designed to construct a "usable and decent" datacube out of the
    mess generated by the emission line fitting function, i.e. out of the many pickle 
    files generated.
    '''
    
    # First, load the original datacube. I need to know how much stuff was fitted.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    # Calculate the outer bin edges for the spectrum
    be = np.append(lams-header1['CD3_3']/2.,lams[-1]+header1['CD3_3']/2.) 
    
    nrows = header1['NAXIS1']
       
    # How many emission lines were fitted ?
    nlines = len(params['elines'].keys()) 
       
    # Very well, what do I want to extract ?
    # 1) A "full emission line spectrum"
    # 2) For each line, a flux map, an intensity map, a velocity map and a dispersion map,
    # and h3 and h4 maps
    # Let's get to it.
    elines_fullspec_cube = data * np.nan # That way, I keep the nan's in the raw data
    # 7 params: Flux, I, v, sigma, h3, h4, SNR
    elines_params_cube = np.zeros((7*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    # And the associated errors !
    elines_params_err = np.zeros((7*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    # Also save the status from mpfit. Could be useful for sorting the good from the bad.
    elines_fit_status = np.zeros((header1['NAXIS2'],header1['NAXIS1']))*np.nan
    
    if params['verbose']:
        print '-> Constructing the datacube for the emission line fitting parameters.'
       
    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cubes [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          fn_list['elines_pickle']+str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            ms = pickle.load(myfile)  
            myfile.close() 
        
            # Get all the parameters in a big array
            ps = [item.params for item in ms]
            errs = [item.perror for item in ms]
            stats = [item.status for item in ms]
            
            # Here, I need to make sure the ps and errs array have a decent shape, even 
            # when the fit failed. Also store the STD that comes out of mpfit (NOT the variance!) 
            ps = [np.zeros_like(ps[0])*np.nan if not(item.status in [1,2,3,4]) else ps[j] 
                                                            for (j,item) in enumerate(ms)]
            
            errs = [np.zeros_like(ps[0])*np.nan if not(item.status in [1,2,3,4]) else 
                                                 errs[j] for (j,item) in enumerate(ms)]
            
            
            # Fill the corresponding datacube
            # First, what are my indices ? I need to keep one slice empty for the SNR of 
            # each line. Get the indices where I will thus store the fit params.
            ind = np.array([np.arange(6)+7*j for j in range(nlines)]).reshape(nlines*6)
            
            elines_params_cube[ind,:,row] = np.array(ps).T
            elines_params_err[ind,:,row] = np.array(errs).T
            elines_fit_status[:,row] = np.array(stats).T
            
            # now, reconstruct the full emission line spectrum
            elines_specs = np.array([brutus_elf.els_spec(lams,p,be=be, 
                                                        method=params['line_profile'],
                                                        inst=params['inst']) for p in ps])
       
            # Fill the corresponding cube 
            elines_fullspec_cube[:,:,row] = elines_specs.T
            
    
    # Now, for each line, the first of these slices is the reference wavelength. 
    # Replace this by the total flux instead ! And make some plots if requested.
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        
        # Calculate the sigma of the fit, in A (incl. instrument dispersion,etc ...)
        # as well as the associated error.
        # WARNING: unlike the data, these are the 1-sigma error coming out of mpfit.
        # Mind the **2 !
        zlams = elines_params_cube[7*k] * (1.+ elines_params_cube[7*k+2] / c )
        zlams_err = elines_params_cube[7*k] * elines_params_err[7*k+2]/c # The 1-std err
        sigma_obs_A = brutus_elf.obs_sigma(elines_params_cube[7*k+3],zlams, 
                                          inst=params['inst'], 
                                          in_errs=[elines_params_err[7*k+3],zlams_err] )
        
        # Compute the line flux
        elines_params_cube[7*k,:,:] = np.sqrt(2*np.pi) * elines_params_cube[7*k+1,:,:] * \
                                      sigma_obs_A[0]                           
        # What about the error ?
        elines_params_err[7*k,:,:] = \
            np.abs(elines_params_cube[7*k,:,:]) * \
            np.sqrt(elines_params_err[7*k+1,:,:]**2/elines_params_cube[7*k+1,:,:]**2 +
                    sigma_obs_A[1]**2/sigma_obs_A[0]**2)           
    
    
    # Now for each line, calculate the SNR
    # For this, I need the continuum cubes
    if fn_list['cont_mix_cube']:
        fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
        hdu = pyfits.open(fn)
        cont_mix_cube = hdu[1].data
        hdu.close()
    else:
        cont_mix_cube = np.zeros_like(data)
        
    # Compute the residual cube
    residuals = data - cont_mix_cube - elines_fullspec_cube
    
    cont_range = params['cont_range']
    lams0 = lams / (params['z_target'] + 1)
    
    # Compute the noise over the region of interest defined by the user
    noise = np.nanstd(residuals[(lams0>=cont_range[0])*(lams0<=cont_range[1]),:,:],axis=0) 
    
    for (k, key) in enumerate(np.sort(params['elines'].keys())):
        elines_params_cube[7*k+6,:,:] = elines_params_cube[7*k+1,:,:]/noise
       
    
    # Export the cube for each emission line parameters as a multi-extension fits file        
    # Do the same for the errors - i.e. params and errors are in two distinct cubes
    for (e,epc) in enumerate([elines_params_cube,elines_params_err]):
        hdu0 = pyfits.PrimaryHDU(None,header0)
    
        hdus = [hdu0]
        # Use the sorted keys, to ensure the same order as the fit parameters
        for (k,key) in enumerate(np.sort(params['elines'].keys())):
            hduk = pyfits.ImageHDU(epc[7*k:7*k+7,:,:])
            # Make sure the WCS coordinates are included as well
            hduk = brutus_tools.hdu_add_wcs(hduk,header1)
            # Also include a brief mention about which version of brutus is being used
            hduk = brutus_tools.hdu_add_brutus(hduk,suffix)
            # Add the line reference wavelength for future references
            hduk.header['BR_REFL'] = (params['elines'][key][0][0], 'reference wavelength')
            hduk.header['BR_CKEY'] = ('F,I,v,sigma,h3,h4,SNR','Content of the cube planes')
            
            hdus.append(hduk)
            
        hdu = pyfits.HDUList(hdus=hdus)
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_elines_'+
                              ['params','perror'][e]+'.fits')
        hdu.writeto(fn_out, clobber=True)
        
        # Add the filename to the dictionary of filenames
        fn_list['elines_'+['params','perror'][e]+'_cube'] = suffix+'_'+params['target']+\
                                                            '_elines_'+\
                                                            ['params','perror'][e]+'.fits'
                   
    # Very well, now let's also create a fits file to save the full emission line spectrum
    # as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(elines_fullspec_cube)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
    hdu1 = brutus_tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_elines_fullspec.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['elines_spec_cube'] = suffix+'_'+params['target']+'_elines_fullspec.fits'
    
    # And finally, the plot with the fit status for each spaxel, to know the bad from the 
    # not-so-bad,
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(elines_fit_status)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_elines_mpfit_status.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['elines_fit_status'] = suffix+'_'+params['target']+'_elines_mpfit_status.fits'
    
    print ' ' 
    return fn_list
# ----------------------------------------------------------------------------------------

def run_plot_elines_cube(fn_list, params, suffix=None, vrange=None, 
                         sigrange=None):   
    '''Creates some plots for the emission lines: F, I, v and sigma.
        
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        vrange: list of int [default: None]
                If set, the range of the colorbar for the velocity plot.
        sigrangre: list of int [default: None]        
                If set, the range of the colorbar for the velocity dispersion plot.
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames.
    '''    
    if params['verbose']:
        print '-> Making some nifty plots from the emission line fitting output.'
       
    fn = os.path.join(params['prod_loc'], fn_list['elines_params_cube'])
    
    # Open the file
    hdu = pyfits.open(fn)
    header0 = hdu[0].header

    # Create a temporary FITS file
    fn_tmp = os.path.join(params['tmp_loc'],suffix+'_tmp.fits')

    # Do this for each emission line fitted
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
    
        # Because of the way aplpy works, I need to stored each "image" in its own fits
        # file. I don't want to keep them, so let's just use a temporary one, get the plot
        # done, and remove it. Not ideal, but who cares ?
        
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
        
        # Make single pretty plots for Flux, Intensity, velocity and velocity dispersion                                        
        for (t,typ) in enumerate(['F','I','v','sigma','h3', 'h4', 'SNR']):
            
            # Take into account the SNR from the user
            if typ != 'SNR':
                this_data[t][this_data[-1] <params['elines'][key][0][3]] = np.nan
            
            # Create a dedicated HDU    
            tmphdu = pyfits.PrimaryHDU(this_data[t])
            # Add the WCS information
            tmphdu = brutus_tools.hdu_add_wcs(tmphdu,this_header)
            tmphdu.writeto(fn_tmp, clobber=True)
            
            # Now plot it 
            this_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                                                        '_eline-'+key+'-'+
                                                        np.str(this_header['BR_REFL'])+
                                                        '_'+
                                                        typ+'.pdf') 
            
            # For the velocity fields, set the vmin and vmax
            if t == 0:
                my_vmin = None
                my_vmax = None
                my_cmap = 'alligator'
                my_stretch = 'arcsinh'
                my_label = r'F [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [125,250,500,1000,2000,4000]
            elif t == 1:
                my_vmin = None
                my_vmax = None
                my_cmap = 'alligator'
                my_stretch = 'arcsinh'
                my_label = r'I [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [50,100,200,400,800,1600]
            elif t == 2 and vrange:
                my_vmin = vrange[0]
                my_vmax = vrange[1]
                my_cmap = 'alligator'
                my_stretch = 'linear'
                my_label = r'$v$ [km s$^{-1}$]'
                my_cbticks = None
            elif t ==3 and sigrange:
                my_vmin = sigrange[0]
                my_vmax = sigrange[1]
                my_cmap = 'alligator'
                my_stretch = 'linear'
                my_label = r'$\sigma_{tot}$ [km s$^{-1}$]'
                my_cbticks = None
            else:
                my_vmin = None
                my_vmax = None
                my_cmap = None
                my_stretch = 'linear'
                my_label = ''
                my_cbticks = None
                                                            
            brutus_plots.make_2Dplot(fn_tmp,ext=0, ofn=this_ofn, contours=False, 
                                    vmin=my_vmin, vmax = my_vmax, cmap=my_cmap,
                                    stretch = my_stretch, cblabel=my_label, 
                                    cbticks = my_cbticks,
                                    scalebar = params['scalebar'],)                                   
    
    # Delete the temporary fits file
    os.remove(fn_tmp)
    # Don't forget to close the initial hdu ...
    hdu.close()
    
    # And also create a plot of the fit status, to see if anything weird happened
    fn = os.path.join(params['prod_loc'], fn_list['elines_fit_status'])
    ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                                                        '_eline_mpfit_status.pdf') 
    brutus_plots.make_2Dplot(fn,ext=1, ofn=ofn, contours=False, vmin=-16, vmax = 8,
                           cmap='alligator',stretch='linear', 
                           cbticks=[-16,0,1,2,3,4,5,6,7,8], cblabel='mpfit status',
                           scalebar = None)
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_inspect_fit(fn_list,params, suffix=None, irange=[None,None], vrange=[None,None]):   
    '''Setup the interactive inspection of the fit results.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        irange: list of int [default: [None,None]]
                The range of the colorbar for the intensity plot.
        vrange: list of int [default: [None,None]]    
                The range of the colorbar for the velocity dispersion plot.

    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames.
    '''
    
    if params['verbose']:
        print '-> Starting the interactive inspection of the fitting.'
    
    # Load the different files I need  
    # Raw data:    
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
    # The Lowess continuum fit
    if not(fn_list['lowess_cube'] is None):
        fn = os.path.join(params['prod_loc'],fn_list['lowess_cube'])
        hdu = pyfits.open(fn)
        lowess = hdu[1].data
        hdu.close()  
    else:
        lowess = np.zeros_like(data)*np.nan
        
    # The ppxf continuum fit
    if not(fn_list['ppxf_cube'] is None):
        fn = os.path.join(params['prod_loc'],fn_list['ppxf_cube'])
        hdu = pyfits.open(fn)
        ppxf = hdu[1].data
        hdu.close()  
    else:
        ppxf = np.zeros_like(data)*np.nan
    
    # The continuum mix datacube
    if not(fn_list['cont_mix_cube'] is None):
        fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
        hdu = pyfits.open(fn)
        cont_mix = hdu[1].data
        hdu.close()  
    else:
        cont_mix = np.zeros_like(data)*np.nan
    
       
    # The elines fit
    if not(fn_list['elines_spec_cube'] is None):
        fn = os.path.join(params['prod_loc'],fn_list['elines_spec_cube'])
        hdu = pyfits.open(fn)
        elines = hdu[1].data
        hdu.close() 
    else:
        elines = np.zeros_like(data)*np.nan 
    
    # Also open the map and vmap to be displayed. 
    # TODO: give user the choice of image on the right
    if not(fn_list['elines_params_cube'] is None):
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
        hdu = pyfits.open(fn)
        map = hdu[1].data[0]
        vmap = hdu[1].data[2]
        hdu.close()
    else:
        map = np.zeros_like(data[0])*np.nan
        vmap = np.zeros_like(data[0])*np.nan    
    
    # A filename if the user decides to save anything.
    my_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+'_fit_inspection_')

    # Launch the interactive plot
    brutus_plots.inspect_spaxels(lams,data,lowess,ppxf,cont_mix, elines, map, vmap, 
                                 irange, vrange,
                                 ofn = my_ofn)    
        
    return fn_list
# ----------------------------------------------------------------------------------------

def run_find_structures(fn_list,params, suffix=None, interactive_mode=True, 
                        automatic_mode=True):   
    ''' 
    This function is designed to identify structures (e.g. HII regions) in the data from
    a 2D image (i.e. an line intensity map), and save them to a pickle file. When
    interactive_mode=True, the user can manually refine the selection. Set
    automatic_mode=False to skip the automatic detection.
    '''
    
    if params['verbose']:
        print '-> Starting the semi-automated procedure for structure identification.'
    
    # Where am I going to save the apertures information ?
    fn_ap = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_ap_list.pkl')
    
    # Add the filename to the dictionary of filenames
    fn_list['ap_list'] = suffix+'_'+params['target']+'_ap_list.pkl'
    
    # Do we want to build apertures based on multiple maps ? Loop, one after the other.
    for key in params['ap_map_lines']:
        print '   Loading eline %s' %key 
        # Is the file already there ? Do I want to overwrite it ? And/or use its content ?
        if os.path.isfile(fn_ap):
            print '    '
            print '   Existing aperture list (%s)' % fn_ap.split('/')[-1]
            print '   Type [a] to load this aperture list and edit it manually'
            print '        [b] to start from scratch, (file will be overwritten!)'
            print '        [c] to do nothing, and continue to the next step'
            while True:
                letter = raw_input()
            
                if letter in ['a','b','c']:
                    break
                else:
                    print '        [%s] unrecognized. Try again.' % letter
        else:
            letter = 'b'
            
        # Open the file, load the list of apertures
        if letter == 'a':
            f = open(fn_ap, 'r')
            start_aps = pickle.load(f)
            f.close()
        elif letter =='b':
            start_aps = None
        elif letter =='c':
            continue
    
        # Now, open the elines param datacube, and extract the Flux map I want to detect
        # stuctures from.

        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
        hdu = pyfits.open(fn)
        header0 = hdu[0].header
        plane = np.sort(params['elines'].keys()).tolist().index(key)
        data = hdu[plane+1].data[0]
        hdu.close()
        
        # Launch the aperture finding routine
        apertures = brutus_plots.build_ap_list(data, start_aps = start_aps, 
                                              radius = params['ap_radius'],
                                              automatic_mode = automatic_mode,
                                              interactive_mode = interactive_mode,
                                              lam = params['elines'][key][0][0],
                                              save_plot = os.path.join(params['plot_loc'],
                                                                       suffix+'_'+
                                                                       params['target']+
                                                                       '_ap_list_'),
                                             )    
        
        # Only if the user wants to save the apertures, do it
        if apertures:
            # Save the results for later use
            f = open(fn_ap, 'w')
            pickle.dump(apertures,f)
            f.close()
        
    return fn_list
# ----------------------------------------------------------------------------------------

def run_make_ap_cube(fn_list, params, suffix=None, do_plot=True):   
    ''' 
    This function is designed to make a cube from a series of apertures (x,y,rs).
    For compativilty with spaxels-by-spaxels analysis codes (incl.brutus), make the cube
    the same size as the original, and repleace each spectra in a given aperture by the
    total aperture spectra. Spaxels outside any apertures are nan's. 
    Assigned spaxels to one aperture only, in order of decreasing flux peak. This makes 
    the data redondant, but will allow for a rapid and direct processing of the resulting
    cube by brutus.
    '''        
    
    if params['verbose']:
        print '-> Constructing the cube with the integrated '+\
              'aperture spectra.'
    
    # Very well, where is the aperture file ?
    fn_ap = os.path.join(params['prod_loc'],fn_list['ap_list'])
    f = open(fn_ap, 'r')
    start_aps = pickle.load(f)
    f.close()
    
    xs,ys,rs = zip(*start_aps)
    
    # I will also need the Flux map from the strongest line - use that set for the 
    # velocity reference of the line fitting
    for key in params['elines'].keys():
        if params['elines'][key][0][0] == params['ref_dv_line']:
            ref_key = key
    
    # Very well, now load the corresponding flux map
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    hdu = pyfits.open(fn)
    fheader0 = hdu[0].header
    plane = np.sort(params['elines'].keys()).tolist().index(ref_key)
    flux = hdu[plane+1].data[0]
    fheader1 = hdu[plane+1].header
    hdu.close()
    
    # I also need to load the raw data cube. Here I ONLY use the raw cube - and NOT the
    # one possibly deredened for galactic extinction.
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Get all the peak intensities associated with each aperture. Needed for sorting them.
    fs = flux[ys,xs]
    
    # Sort them in decreasing order of peak intensity
    sort_index = np.argsort(fs)[::-1]
    
    # Now, construct a map where each pixel contains the number of the region it belongs 
    # to. Starting from 0.
    ap_map = np.zeros_like(flux) * np.nan
    (ny,nx) = np.shape(flux)
    indices = np.mgrid[0:ny,0:nx]
    
    # Also construct an aperture spectra cube
    ap_spec_cube = np.zeros_like(data) * np.nan
    ap_spec_err = np.zeros_like(data) * np.nan
    
    # Loop through each ap. Hopefully not too slow ...
    for (i,ind) in enumerate(sort_index):
        progress = 100. * (i+1.)/len(sort_index)
        sys.stdout.write('\r   Dealing with aperture %i [%5.1f%s]' % 
                         (i,progress,'%'))
        sys.stdout.flush()
        x = xs[ind]
        y = ys[ind]
        r = rs[ind]
        # Find all spaxels with the ap radius 
        # Don't do anything fancy, just measure the distance to each spaxel center.
        in_ap = (indices[1]-x)**2+(indices[0]-y)**2 <= r**2
        
        # Assign each spaxel (not yet assigned to another brighter feature) the 
        # corresponding ap number.
        ap_map[in_ap * np.isnan(ap_map)] = i

        #For this aperture, sum all spaxels into master aperture spectra, and fill the 
        # cube. Avoid the nans (e.g. mosaic edges, etc ...)

        spec = np.nansum(data[:,indices[0][ap_map==i],indices[1][ap_map==i]],axis=1)
        err = np.nansum(error[:,indices[0][ap_map==i],indices[1][ap_map==i]],axis=1)

        # Loop through each spaxel in the aperture. There MUST be a smarter way, but I
        # can't figure it out. Should not be too costly time-wise anyway ...
        for k in range(len(indices[1][ap_map==i])):
            xi = indices[1][ap_map==i][k]
            yi = indices[0][ap_map==i][k]
            ap_spec_cube[:,yi,xi] = spec
            ap_spec_err[:,yi,xi]  = err
        
    # All done. Save the aperture map to a fits file.
    hdu0 = pyfits.PrimaryHDU(None,fheader0)
    hdu1 = pyfits.ImageHDU(ap_map)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,fheader1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_ap_map.fits')
    hdu.writeto(fn_out, clobber=True)    
    
    # Add this filename to the dictionary of filenames
    fn_list['ap_map'] = suffix+'_'+params['target']+'_ap_map.fits'
    
    # Make a plot of the apertures ?
    if do_plot:
        this_ofn = os.path.join(params['plot_loc'],
                          suffix+'_'+params['target']+'_ap_map.pdf')
        brutus_plots.make_2Dplot(fn_out,ext=1, ofn=this_ofn, contours=False, 
                                    vmin=0, vmax = len(xs), cmap='alligator',
                                    stretch = 'linear', 
                                    cblabel=r'Aperture idendification number', 
                                    cbticks = None,
                                    scalebar = None,
                                )  
    
    # And also save the aperture spectral cube to a fits file
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(ap_spec_cube)
    hdu2 = pyfits.ImageHDU(ap_spec_err)
    # Make sure the WCS coordinates are included as well
    for hdu in [hdu1,hdu2]:
        hdu = brutus_tools.hdu_add_wcs(hdu,header1)
        hdu = brutus_tools.hdu_add_lams(hdu,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu = brutus_tools.hdu_add_brutus(hdu,suffix)
        
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1,hdu2])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_ap_spec_cube.fits')
    hdu.writeto(fn_out, clobber=True)   
      
    # Add this filename to the dictionary of filenames
    fn_list['ap_spec_cube'] = suffix+'_'+params['target']+'_ap_spec_cube.fits'
                                     
    print ' '
    
    return fn_list           
# ----------------------------------------------------------------------------------------

def run_extragal_dered(fn_list, params, suffix=None, do_plot=True):   
    '''Corrects the line fluxes for the extragalactic reddening using Ha and Hb.
    
    This function returns a corrected set of line fluxes, as well as the associated Av 
    map.
     
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make a plot of the Av map or not.
        
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        Should I add some info about each curve here ? 
    '''  
    
    if params['verbose']:
        print '-> Correcting for the extragalactic attenuation.'
    
    # Open the raw data cube
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    nlines = len(params['elines'].keys())

    # Import the Halpha and Hbeta maps, and mask any spaxel not within the proper SNR range
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    hdu = pyfits.open(fn)
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        if params['elines'][key][0][0] == ha:
            ha_map = hdu[k+1].data[0]
            ha_map[hdu[k+1].data[-1] < params['elines'][key][0][3]] = np.nan
        elif params['elines'][key][0][0] == hb:
            hb_map = hdu[k+1].data[0]
            hb_map[hdu[k+1].data[-1] < params['elines'][key][0][3]] = np.nan
    hdu.close()


    hahb = ha_map / hb_map
    # Construct the Av map, just because I can
    av = brutus_red.hahb_to_av(hahb, params['hahb_0'], curve = params['egal_curve'], 
                              rv = params['egal_rv'], rva = params['egal_rva'])

    # Now, for each emission line fitted, let's correct the flux:
    # A storage structure for the emission lines. Keep the un-reddened flux as well.
    drelines_params_cube = np.zeros((8*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    drelines_perror_cube = np.zeros((8*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan

    # Load the current line parameters, loop through, and save to a new - bigger - array.
    # Also take care of the error
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    fn_err = os.path.join(params['prod_loc'],fn_list['elines_perror_cube'])
    hdu = pyfits.open(fn)
    hdu_err = pyfits.open(fn_err)
    fheader0 = hdu[0].header
    fheader0_err = hdu_err[0].header
    
    # Loop through each emission line
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        drelines_params_cube[8*k+1:8*(k+1)] = hdu[k+1].data
        drelines_perror_cube[8*k+1:8*(k+1)] = hdu_err[k+1].data
        # And add the de-reddened line flux in the first layer.
        # Note: the reddening correction is based on the de-redshifted line wavelength,
        # because this happens in the rest-frame of the target !
        this_lam = params['elines'][key][0][0]
        drelines_params_cube[8*k] = hdu[k+1].data[0] * \
                                     brutus_red.extragalactic_red(this_lam, hahb, 
                                                                 params['hahb_0'],
                                                                 curve = params['egal_curve'], 
                                                                 rv = params['egal_rv'],
                                                                 rva = params['egal_rva'])
        # Assume the reddening correction is error free ... sigh...
        drelines_perror_cube[8*k] = hdu_err[k+1].data[0] * \
                                     (brutus_red.extragalactic_red(this_lam, hahb, 
                                                                 params['hahb_0'],
                                                                 curve = params['egal_curve'], 
                                                                 rv = params['egal_rv'],
                                                                 rva = params['egal_rva']))**2                                                        
    fheader1 = hdu[1].header
    fheader1_err = hdu_err[1].header
    hdu.close()
    hdu_err.close()
    
    
    for (e,epc) in enumerate([drelines_params_cube, drelines_perror_cube]):
        # Ok, now save the data 
        hdu0 = pyfits.PrimaryHDU(None,fheader0)
    
        hdus = [hdu0]
        # Use the sorted keys, to ensure the same order as the fit parameters
        for (k,key) in enumerate(np.sort(params['elines'].keys())):
            hduk = pyfits.ImageHDU(epc[8*k:8*(k+1),:,:])
            # Make sure the WCS coordinates are included as well
            hduk = brutus_tools.hdu_add_wcs(hduk,header1)
            # Also include a brief mention about which version of brutus is being used
            hduk = brutus_tools.hdu_add_brutus(hduk,suffix)
            # Add the line reference wavelength for future references
            hduk.header['BR_REFL'] = (params['elines'][key][0][0], 'reference wavelength')
            hduk.header['BR_CKEY'] = ('dredF,F,I,v,sigma,h3,h4','Content of the cube planes')
            
            hdus.append(hduk)
            
        hdu = pyfits.HDUList(hdus=hdus)
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_elines_'+
                              ['params','perror'][e]+'_dered.fits')
        hdu.writeto(fn_out, clobber=True)
        
        # Add the filename to the dictionary of filenames
        fn_list['dered_elines_'+['params','perror'][e]] = suffix+'_'+\
                                                                  params['target']+\
                                                                  '_elines_'+\
                                                                  ['params','perror'][e]+\
                                                                  '_dered.fits'
                   
    # Very well, now let's also create a fits file to save the Av map
    # as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(av)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
    # Also include a brief mention about which version of brutus is being used
    hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_Av.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['Av_map'] = suffix+'_'+params['target']+'_Av.fits'
    
    # If requested, save a plot for the Av map.
    if do_plot:
        
        this_ofn = os.path.join(params['plot_loc'],
                          suffix+'_'+params['target']+'_Av_map.pdf')
        brutus_plots.make_2Dplot(fn_out,ext=1, ofn=this_ofn, contours=False, 
                                    vmin=0, vmax = 2, cmap='alligator',
                                    stretch = 'linear', 
                                    cblabel=r'A$_V$ [mag]', 
                                    cbticks = None,
                                    scalebar = params['scalebar'],
                                )
                                    
    return fn_list
# ---------------------------------------------------------------------------------------- 
  
def run_fit_kinematic_pa(fn_list, params, suffix = None, do_plot = True, 
                         vrange = [None,None]):
    '''
    This function computes the kinematic position angle for the stars or the gas, 
    via the fit_kinematic_pa() routine from M.Cappellari (not included with brutus). 
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make a plot of the kinematic map or not.
        vrange: list of int [default: [None,None]]    
                The range of the colorbar for the velocity plot.
                 
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        pyqz is NOT included with brutus. A separate installation is required.
        See http://fpavogt.github.io/pyqz/installation.html
    '''
    
    if params['verbose']:
        print '-> Fitting the kinematic position angle.'
    
    all_kin_pas = {}
    
    # How many kinematic maps do I want to fit ?
    for kin in params['fit_kin']:
        if params['verbose']:
            print '   [%s] ...' % kin
            
        # First, I need to load the kinematic maps and the associated errors
        if kin == 'ppxf': 
        
            fn = os.path.join(params['prod_loc'],fn_list['ppxf_sol_map']) 
            hdu = pyfits.open(fn)
            header1 = hdu[1].header
            vel = hdu[1].data[0]
            dvel = None # No error on the velocity coming from ppxf. Will be set to 10km/s 
                        # by fit_kinematic_pa()
            hdu.close()
            
        elif kin in params['elines'].keys():
            
            ind = np.sort(params['elines'].keys()).tolist().index(kin)
            
            fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
            hdu = pyfits.open(fn)
            header1 = hdu[ind+1].header
            vel = hdu[ind+1].data[2]
            snr = hdu[ind+1].data[-1]
            vel[snr<params['elines'][kin][0][3]] = np.nan # Take into account the SNR
            hdu.close()
            
            fn = os.path.join(params['prod_loc'],fn_list['elines_perror_cube'])
            hdu = pyfits.open(fn)
            dvel = hdu[ind+1].data[2]
            hdu.close()
        
        # Correct for the redshift
        vel -= params['z_target'] * c # Correct for the velocity
                
        # Very well, I now have my velocity map. How big is it ?
        (ny,nx) = np.shape(vel)
        
        # Build the array of indices
        ys,xs = np.mgrid[0:ny:1,0:nx:1]
        
        # Now, I need to recenter this based on the kinematical center for the object
        ys = ys.astype(float) - params['fit_kin_center'][1]
        xs = xs.astype(float) - params['fit_kin_center'][1]
        
        # Handle the presence of NaNs ...
        xs = xs[vel==vel]
        ys = ys[vel==vel]
        vel_nonan = vel[vel==vel]
        if not(dvel is None):
            dvel = dvel[vel==vel]
    
        # Reschape all of it ...    
        xs = xs.reshape(np.size(xs))
        ys = ys.reshape(np.size(ys))
        vel_nonan = vel_nonan.reshape(np.size(vel_nonan))
        if not(dvel is None):
            dvel = dvel.reshape(np.size(dvel))
            
        # Very well, let's launch the processing
        all_kin_pas[kin] = fkp.fit_kinematic_pa(xs, ys, vel_nonan, dvel = dvel,
                                                debug = False, 
                                                nsteps = params['fit_kin_nsteps'], 
                                                plot = False)
    
        # Ok, what do I want to do with this ? A plot, and save it for later ?
        if do_plot:
            
            # Create a temporary FITS file
            fn_tmp = os.path.join(params['tmp_loc'],suffix+'_tmp.fits')
            
            # Let's create a temporary fits file with the given velocity map  
            tmphdu = pyfits.PrimaryHDU(vel + params['z_target'] * c)
            # Add the WCS information
            tmphdu = brutus_tools.hdu_add_wcs(tmphdu,header1)
            tmphdu.writeto(fn_tmp, clobber=True)
            
        
            # Alright, let's make a velocity plot.
            if kin in params['elines'].keys():
                ofn = suffix+'_'+params['target']+'_fit_kinematic_pa_%s.pdf' % \
                                                              params['elines'][kin][0][0]
            else:
                ofn = suffix+'_'+params['target']+'_fit_kinematic_pa_%s.pdf' % kin
            
            ofn = os.path.join(params['plot_loc'], ofn)
            
            brutus_plots.make_2Dvelplot(fn_tmp, # path to the data (complete!)
                                        ext = 0, # Which extension am I looking for ?
                                        ofn = ofn, # Savefig filename
                                        contours = False, # Draw any contours 
                                        stretch = 'linear',
                                        vmin = vrange[0], 
                                        vmax = vrange[1],
                                        cmap = 'alligator',
                                        cblabel = r'v [km s$^{-1}$]',
                                        cbticks = None,
                                        scalebar = params['scalebar'],
                                        pa = all_kin_pas[kin][:2], 
                                        center = params['fit_kin_center'],
                                        )
            # Clean up the temporary file
            os.remove(fn_tmp)
            
    # Add the filename to the dictionary of filenames
    fn_list['fit_kin_pa'] = suffix+'_'+params['target']+'_fit_kinematic_pa.pkl'
    
    # And save these avlues in a pickle file, for future access if needed.
    f = open(os.path.join(params['prod_loc'],fn_list['fit_kin_pa']), 'w')
    pickle.dump(all_kin_pas, f)
    f.close()

    return fn_list
# ---------------------------------------------------------------------------------------- 
  
def run_get_ne(fn_list, params, suffix = None, do_plot = True, ratio_range=[None, None]):
    '''
    This function computes the kinematic position angle for the stars or the gas, 
    via the fit_kinematic_pa() routine from M.Cappellari (not included with brutus). 
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                Whether to make a plot of the electron density map or not.
        ratio_range: list [default[None,None]
                The line ratio range, for the plot.
                 
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        pyqz is NOT included with brutus. A separate installation is required.
        See http://fpavogt.github.io/pyqz/installation.html
    '''
    
    if params['verbose']:
        print '-> Constructing the [SII] line ratio.'
    
    
    # Alright, let's load the emission line fluxes
    if params['pyqz_use_egal_dered']:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
        fn_e = os.path.join(params['prod_loc'],fn_list['dered_elines_perror'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
        fn_e = os.path.join(params['prod_loc'],fn_list['elines_perror_cube'])

    
    s2h_key = params['elines']
    
    # Open the files
    hdu = pyfits.open(fn)
    header0 = hdu[0].header

    hdu_e = pyfits.open(fn_e)
    header0_e = hdu_e[0].header
    
    
    # Now loop through these, until I find the lines I need
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        
        if params['elines'][key][0][0] == s2h:
            header_s2h = hdu[k+1].header
            s2h_map = hdu[k+1].data[0]
            s2h_err = hdu_e[k+1].data[0]
        elif params['elines'][key][0][0] == s2l:
            header_s2l = hdu[k+1].header
            s2l_map = hdu[k+1].data[0]
            s2l_err = hdu_e[k+1].data[0]
    
    # Very well, let's compute the ratio map and error:
    ratio_map = s2h_map/s2l_map
    ratio_err = np.abs(ratio_map) * np.sqrt((s2h_err/s2h_map)**2 + (s2l_err/s2l_map)**2)
    
    # Let's export this to a dedicated fits file.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(ratio_map)
    hdu2 = pyfits.ImageHDU(ratio_err)
    for hdu in [hdu1,hdu2]:
        # Make sure the WCS coordinates are included as well
        hdu = brutus_tools.hdu_add_wcs(hdu, header_s2h)
        # Also include a brief mention about which version of brutus is being used
        hdu = brutus_tools.hdu_add_brutus(hdu, suffix)
    
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1,hdu2])
    fn_out = os.path.join(params['prod_loc'], suffix+'_'+params['target']+'_[SII]-ratio.fits')
    hdu.writeto(fn_out, clobber=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['ne_map'] = suffix+'_'+params['target']+'_[SII]-ratio.fits'
    
    # And now, let's plot it.
    if do_plot:
         
            ofn = suffix+'_'+params['target']+'_[SII]-ratio.pdf'
            
            ofn = os.path.join(params['plot_loc'], ofn)
            
            brutus_plots.make_2Dplot(fn_out, # path to the data (complete!)
                                        ext = 1, # Which extension am I looking for ?
                                        ofn = ofn, # Savefig filename
                                        contours = False, # Draw any contours 
                                        stretch = 'linear',
                                        vmin = ratio_range[0], 
                                        vmax = ratio_range[1],
                                        cmap = 'alligator',
                                        cblabel = r'[S\,\textsc{\smaller{II}}] $\lambda$6731 / [S\,\textsc{\smaller{II}}] $\lambda$6716 ',
                                        cbticks = None,
                                        scalebar = params['scalebar'],
                                        ) 
        
    
    return fn_list

# ---------------------------------------------------------------------------------------- 
  
def run_get_QZ(fn_list, params, suffix = None, start_row = 0, end_row = None):
    '''
    This function computes the ionization parameter and oxygen abundance of HII regions
    using the pyqz module: http://fpavogt.github.io/pyqz/ 
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        start_row: int [default: 0]
                   The starting row from which to process the data
        end_row: int [default: None]
                 The row at which to finish the processing. Set to None for the max row.
                 
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        pyqz is NOT included with brutus. A separate installation is required.
        See http://fpavogt.github.io/pyqz/installation.html
    '''       
    
    if params['verbose']:
        print '-> Deriving Log(Q) and Tot(O)+12 using pyqz.'
    
    # Very well, now I need to open the line fluxes.
    if params['pyqz_use_egal_dered']:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
        fn_e = os.path.join(params['prod_loc'],fn_list['dered_elines_perror'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
        fn_e = os.path.join(params['prod_loc'],fn_list['elines_perror_cube'])

    # Alright, now I need to extract the line fluxes and errors
    # Open the files
    hdu = pyfits.open(fn)
    header0 = hdu[0].header

    hdu_e = pyfits.open(fn_e)
    header0_e = hdu_e[0].header

    # All elines definition are inside elines_pyqz inside brutus_metadata. 
    # For each of them, construct a flux map.

    elines_fluxes = {}
    elines_errors = {}

    # Look at each emission line fitted. If it is part of a "pyqz emission line
    for (k,akey) in enumerate(np.sort(params['elines'].keys())):
    
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
        this_error = hdu_e[k+1].data
                               
        lam = this_header['BR_REFL']
        
        for key in elines_pyqz.keys(): # Search all possible known pyqz lines
            if lam in elines_pyqz[key]: # Wavelength match ? I found a line.
            
                flux = this_data[0]
                err = this_error[0]
                
                # Take into account the SNR
                flux[this_data[-1]<params['elines'][akey][0][3]] = np.nan
            
                # Already exists ? - then I found a line with multiple components
                if key in elines_fluxes.keys(): 
                    elines_fluxes[key].append(flux) # append the flux
                    elines_errors[key].append(err)
                else:
                    elines_fluxes[key] = [flux] # Save the data for later
                    elines_errors[key] = [err]
                       
    # Close the hdus
    hdu.close()
    hdu_e.close()
                    
    # Alright, I extracted all my lines fluxes. Am I missing anything ? 
    for key in elines_fluxes.keys():
        if len(elines_fluxes[key]) != len(elines_pyqz[key]):
            sys.exit('ERROR: a line component is missing for %s. Was it fitted ?' % key)

    # Now, proceed to construct the dictionary of line fluxes and errors.
    nlines = len(elines_fluxes.keys())

    nx = this_header['NAXIS1']
    ny = this_header['NAXIS2']

    # Get some info about the cube
    nrows = this_header['NAXIS1']
    if start_row is None:
        start_row = 0
    if end_row is None:
	    end_row = nrows-1 
    
    # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
	# loose the current row. Just live with it.
	# Here, multiprocessing is already implemented inside the pyqz module. No need to do
	# anything else but feed pyqz.get_global_qz() big chunks in a for loop.
    for row in np.linspace(start_row, end_row, end_row-start_row+1):       
	    
        if params['multiprocessing']:
            if type(params['multiprocessing']) == np.int:    
                nproc = params['multiprocessing']
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
        else:
            nproc = 1
        
        # Limit this to 8 for brutus - memory issues ...
        if nproc>8:
            nproc = 8    
	        
        sys.stdout.write('\r   Dealing with spectra in row %2.i, %i at a time ...' % 
                         (row,nproc))
        sys.stdout.flush() 
        
        # The structure to hold the list of line fluxes, for each row
        allmylines = np.zeros((ny,nlines*2))
        allmynames = []
        
        # Very well, let's start filling a big array with all my lines fluxes lined up 
        for (k,key) in enumerate(elines_fluxes.keys()):

             # Very well, sum it all
            fluxes = np.sum(elines_fluxes[key],axis=0)
            errors = np.sqrt(np.sum([item**2 for item in elines_errors[key]],axis=0)) # CHECK THAT THIS DOES WHAT I THINK IT DOES ! 
            # pyqz requires errors as 1-std. Like the output of mpfit.
            # Also get rid of the nan's in the errors ... because pyqz can't deal with them.
            errors[np.isnan(errors)] = 0.0

            allmylines[:,2*k] = fluxes[:,row].reshape(ny)
            allmylines[:,2*k+1] = errors[:,row].reshape(ny)
    
            allmynames.append(key)
            allmynames.append('std'+key)
    
        # I need to keep track of which spaxel went where, so I can make sense of the files
        # generated by pyqz. 
        spaxel_ids = ['%i-%i' % (row,i) for i in range(ny)]    
    
        # Ready to go - let's launch pyqz
        pyqz_out = pyqz.get_global_qz(allmylines, allmynames, params['pyqz_diags'], 
                                      ids = spaxel_ids, 
                                      qzs = ['LogQ', 'Tot[O]+12'],
                                      Pk=5.0, kappa = np.inf, struct='pp', sampling=1, 
                                      error_pdf = 'normal', 
                                      srs=400, flag_level=2.0,  KDE_method='gauss', 
                                      KDE_qz_sampling = 101j, KDE_do_singles = True, 
                                      KDE_pickle_loc = params['pyqz_loc'], 
                                      verbose = False, nproc = nproc)    
 
	    # Here, I need to save these results. Pickle could be fast and temporary,
	    # Until I then re-build the entire cube later on ? Also allow for better
	    # row-by-row flexibility.
        fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_pyqz_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')
        file = open(fn,'w')
        pickle.dump(pyqz_out,file)
        file.close()
   	
   	# Add the generic filename to the dictionary of filenames
   	fn_list['pyqz_pickle'] = suffix+'_'+params['target']+'_pyqz_row_'
   	     
    print ' done !'
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_make_QZ_cube(fn_list, params, suffix=None, do_plot = True):   
    ''' 
    This function is designed to construct a "usable and decent" datacube out of the
    mess generated by the pyqz step.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make some nifty plots or not.
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
    ''' 
    
    # First, load the original datacube. I need to know how much stuff was fitted.
    if fn_list['galdered_cube'] is None:
        hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = pyfits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
        
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close() 
    
    nrows = header1['NAXIS1']
       
    # Very well, what do I want to extract ?
    # A 2-D map for all the variable returned by pyqz. STore this in a multi-extension fits
    
    # Open one of these files, to check how much variables I want to map
    fns = glob.glob(os.path.join(params['tmp_loc'],'*pyqz*.pkl'))
    
    ftmp = open(fns[0],'r')
    content = pickle.load(ftmp)
    ftmp.close()
    
    map_names = content[1] # All the different variable names.
   
    all_maps = np.zeros((len(map_names),header1['NAXIS2'],header1['NAXIS1']))*np.nan
    
    if params['verbose']:
        print '-> Constructing the datacube for the pyqz output.'
       
    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cube [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          fn_list['pyqz_pickle']+str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            ms = pickle.load(myfile)  
            myfile.close() 
        
            # Ok, one line to store it all back !
            all_maps[:,:,row] = ms[0].T
    
    # Get ready with the extension file
    hdu0 = pyfits.PrimaryHDU(None,header0)
    
    hdus = [hdu0]
    # Store each map in a dedicated HDU
    for (k,key) in enumerate(map_names):
        hduk = pyfits.ImageHDU(all_maps[k,:,:])
        # Make sure the WCS coordinates are included as well
        hduk = brutus_tools.hdu_add_wcs(hduk,header1)
        # Also include a brief mention about which version of brutus is being used
        hduk = brutus_tools.hdu_add_brutus(hduk,suffix)
        # Add the line reference wavelength for future references
        hduk.header['BR_PYQZ'] = (key, 'pyqz variable')
            
        hdus.append(hduk)
            
    hdu = pyfits.HDUList(hdus=hdus)
    fn_out = os.path.join(params['prod_loc'], 
                          suffix+'_'+params['target']+'_pyqz_QZs.fits')
    hdu.writeto(fn_out, clobber=True)
        
    # Add the filename to the dictionary of filenames
    fn_list['pyqz_QZs'] = suffix+'_'+params['target']+'_pyqz_QZs.fits'
    
    print ' '
    
    # Make some plots, if requested:
    for (m,map) in enumerate(map_names):
        
        
        ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+'_pyqz_'+
                         map.replace(';','_vs_').replace('|','_').replace('/','-')+'.pdf')
        
        mymap = map.replace('<','$<$').replace('>','$>$').replace('{','$\{$')
        mymap = mymap.replace('}','$\}$').replace('Log','$\log$').replace('|','  ')
        mymap = mymap.replace('_','\_')
        mylabel = r'' + mymap
        
        brutus_plots.make_2Dplot(fn_out, # path to the data (complete!)
                    ext = m+1, # Which extension am I looking for ?
                    ofn = ofn, # Savefig filename
                    contours = False, # Draw any contours 
                    stretch = 'linear',
                    vmin = None, 
                    vmax = None,
                    cmap = 'alligator',
                    cblabel = mylabel,
                    cbticks = None,
                    scalebar = params['scalebar'],
                )
    
    return fn_list    
# ----------------------------------------------------------------------------------------

def run_plot_elines_RGB(fn_list, params, suffix=None, 
                        mixes = [['[NII]','Ha','[OIII]'],], 
                        stretches = ['log',],
                        stretch_plims = [[10.,99.5,10.,99.5,10.,99.5],],
                        stretch_vlims = [[None,None,None,None,None,None],],
                        use_egal_dered = False,
                       ):   
    ''' 
    This function is designed to make some RGB images from emission lines.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        mixes: list of list [default: [['[NII]','Ha','[OIII]']]]
             A list of line triplets indicatin the R, G & B channels.
        stretches: list of string [default: ['log']]
                   The stretches to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
        stretch_plims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                       The limiting percentiles for the plot, as 
                       [pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b]
        stretch_vlims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                       The limtiing values for the plot (superseeds stretch_plims), as
                       [vmin_r, vmax_r, vmin_g, vmax_g, vmin_b, vmax_b]
        use_egal_dered: bool [default: False]
                        If available, whether to use the line fluxes corrected for 
                        extragalactic attentuation.
             
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
    ''' 
    
    # First, I need to generate 3 individual fits files with the line fluxes
    if params['verbose']:
        print '-> Creating some RGB emission line images.'
    
    # Very well, now I need to open the line fluxes.
    if use_egal_dered:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])

    # Alright, now I need to extract the line fluxes
    # Open the file
    hdu = pyfits.open(fn)
    header0 = hdu[0].header
    header1 = hdu[1].header

    # All elines definition are inside elines_pyqz inside brutus_metadata. 
    # For the ones the user requested, construct a flux map.
    elines_fluxes = {}
    for (m, mix) in enumerate(mixes):
        elines_fluxes[m] = {0:None, 1:None, 2:None}    

    # Look at each emission line fitted. Do I need it ?
    for (k,akey) in enumerate(np.sort(params['elines'].keys())):
    
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
                               
        lam = this_header['BR_REFL']
        
        for (m,mix) in enumerate(mixes):
            
            for (cc,col) in enumerate(mix): # Search all possible lines
                if lam in elines_pyqz[col]: # Wavelength match ? I found a line.
            
                    flux = this_data[0]
                
                    # Take into account the SNR
                    flux[this_data[-1]<params['elines'][akey][0][3]] = np.nan
            
                    # Already exists ? - then I found a line with multiple components
                    if not(elines_fluxes[m][cc] is None): 
                        elines_fluxes[m][cc] += flux # append the flux
                    else:
                        elines_fluxes[m][cc] = flux # Save the data for later
                       
    # Close the hdus
    hdu.close()
    
    # Very well, I now have all my line fluxes. Let's deal with each RGB set one at a time.
    for (m,mix) in enumerate(mixes):

        fns = []
        
        ofn = os.path.join(params['plot_loc'],
                 suffix+'_'+params['target']+'_RGB_%s-%s-%s.pdf' % (mix[0],mix[1],mix[2]))
                 
        # Now, let' create 3 independant temporary fits file with WCS coordinates
        for (cc,col) in enumerate(mix):
            fn = os.path.join(params['plot_loc'],'RGB_tmp_%i.fits' % cc)
            fns.append(fn)
            
            hdu = pyfits.PrimaryHDU(elines_fluxes[m][cc], header = header1)
            # Add the wcs info
            hdu = brutus_tools.hdu_add_wcs(hdu, header1)
            outfits = pyfits.HDUList([hdu])
            outfits.writeto(fn,clobber=True)
    
        # Great, I am now ready to call the plotting function
        brutus_plots.make_RGBplot(fns, ofn,  stretch = stretches[m],
                                  plims = stretch_plims[m], vlims = stretch_vlims[m],
                                  title = r'RGB: %s %s %s' %  (mix[0],mix[1],mix[2]),
                                  scalebar = params['scalebar'])
    
        # And remember to delete all the temporary files
        for f in fns:
            os.remove(f)    
    
    return fn_list    
# ----------------------------------------------------------------------------------------

def run_plot_flux_ratio(fn_list, params, suffix=None, 
                        ratios = ['[NII]/[OIII]',],
                        vrange = [[None,None],],
                        use_egal_dered = False,
                       ):   
    ''' 
    This function is designed to make some images from emission lines flux ratios.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        ratios: list of string [default: ['[NII]/[OIII]']]
             A list of line ratios to plot.
        vrange: list of lists [default: [[Mone,None]]]
                The plot range.
        use_egal_dered: bool [default: False]
                        If available, whether to use the line fluxes corrected for 
                        extragalactic attentuation.
             
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
    ''' 
    
       # First, I need to generate 3 individual fits files with the line fluxes
    if params['verbose']:
        print '-> Creating some line ratio images.'
    
    # Very well, now I need to open the line fluxes.
    if use_egal_dered:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])

    # Alright, now I need to extract the line fluxes
    # Open the file
    hdu = pyfits.open(fn)
    header0 = hdu[0].header
    header1 = hdu[1].header

    # All elines definition are inside elines_pyqz inside brutus_metadata. 
    # For the ones the user requested, construct a flux map.
    elines_fluxes = {}
    for (m, mix) in enumerate(ratios):
        elines_fluxes[m] = {0:None, 1:None}    

    # Look at each emission line fitted. Do I need it ?
    for (k,akey) in enumerate(np.sort(params['elines'].keys())):
    
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
                               
        lam = this_header['BR_REFL']
        
        for (m,mix) in enumerate(ratios):
            
            for (cc,col) in enumerate(mix.split('/')): # Search all possible lines
                if lam in elines_pyqz[col]: # Wavelength match ? I found a line.
            
                    flux = this_data[0]
                
                    # Take into account the SNR
                    flux[this_data[-1]<params['elines'][akey][0][3]] = np.nan
            
                    # Already exists ? - then I found a line with multiple components
                    if not(elines_fluxes[m][cc] is None): 
                        elines_fluxes[m][cc] += flux # append the flux
                    else:
                        elines_fluxes[m][cc] = flux # Save the data for later
                       
    # Close the hdus
    hdu.close()
    
    # Very well, I now have all my line fluxes. Let's deal with each set one at a time.
    for (m,mix) in enumerate(ratios):
        
        ofn = os.path.join(params['plot_loc'],
                 suffix+'_'+params['target']+'_ratio_%s_vs_%s.pdf' % (mix.split('/')[0],
                                                                      mix.split('/')[1]))
                 
        # Now, let' create a fits file with WCS coordinates
        fn = os.path.join(params['prod_loc'],
                 suffix+'_'+params['target']+'_ratio_%s_vs_%s.fits' % (mix.split('/')[0],
                                                                       mix.split('/')[1]))
            
        hdu0 = pyfits.PrimaryHDU(header = header0)
        hdu1 = pyfits.ImageHDU(np.log10(elines_fluxes[m][0]/elines_fluxes[m][1]))
        # Make sure the WCS coordinates are included as well
        hdu1 = brutus_tools.hdu_add_wcs(hdu1,header1)
        # Also include a brief mention about which version of brutus is being used
        hdu1 = brutus_tools.hdu_add_brutus(hdu1,suffix)
        # Add the line reference wavelength for future references
        hdu1.header['BR_RATIO'] = (mix, 'line ratio')
            
        hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
        hdu.writeto(fn, clobber=True)
        
        # Great, I am now ready to call the plotting function
        brutus_plots.make_2Dplot(fn, ext = 1, ofn = ofn, contours = False, 
                                 stretch = 'linear',
                                 vmin = vrange[m][0], 
                                 vmax = vrange[m][1],
                                 cmap = alligator,
                                 cblabel = r'$\log$ %s' %mix,
                                 cbticks = None,
                                 scalebar = params['scalebar'],
                                )
                                
    
    return fn_list  