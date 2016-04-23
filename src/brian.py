# -*- coding: utf-8 -*-
#
# This file contains BRIAN routines to fit the stellar continuum and the emission lines
# in an IFU data cube (i.e. MUSE).
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
import sys
import os
from astropy.io import fits as pyfits
from functools import partial
import pickle
import multiprocessing

from brian_metadata import __version__

import brian_tools as tools
import brian_cof
import brian_elf
import brian_plots
from brian_metadata import *
   
# ---------------------------------------------------------------------------------------- 
  
def run_snr_maps(params, suffix = None, do_plot = False):
    '''
    This function computes the SNR maps for the continuum and Ha (or other line) for a 
    MUSE datacube. It also creates a map of spaxels with any signal at all.
    The resulting maps are saved to a fits file with full header and WCS coordinates.
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
    # The signal
    cont_s = np.median(data[(lams>=cont_range[0])*(lams<=cont_range[1]),:,:],axis=0) 
    # The noise
    cont_n = np.std(data[(lams>=cont_range[0])*(lams<=cont_range[1]),:,:],axis=0) 
    
    cont_snr = cont_s/cont_n
    # Also make sure this is always > 0
    cont_snr[cont_snr <0] = 0
    
    # I also want to compute the SNR map for the strongest emission line in there.
    # This line is defined via the 'ref_dv_line' in params (i.e. it is the same line used 
    # for initial v guess for the emission line fitting).
    line_s = np.max(data[ (lams>=line_range[0]) * (lams<=line_range[1]), :,:], axis=0)
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
        hdu = tools.hdu_add_wcs(hdu,header1)
        # Also include a brief mention about which version of BRIAN is being used
        hdu.header['BRIAN_V'] = (__version__,'brian version that created this file.')
    # For reference, also include the line/region this maps are based on
    hdu1.header['BRIAN_R'] = (np.str(cont_range), 'spectral range used for continuum SNR')
    hdu2.header['BRIAN_R'] = (np.str(line_range), 'spectral range used for line SNR') 
        
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1,hdu2,hdu3])
    fn_out = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_snr.fits')
    hdu.writeto(fn_out, clobber=True)
                
    if do_plot:
        # Alright, let's take out the big guns ...            
        brian_plots.make_2Dplot(fn_out,1, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_cont_snr.pdf'),
                            contours = [3,5,10,20], vmin=0,vmax=30,
                            )  
        brian_plots.make_2Dplot(fn_out,2, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_line_snr.pdf'),
                            contours = [3,5,10,20], vmin=0,vmax=30,
                            )    
        brian_plots.make_2Dplot(fn_out,3, 
                            os.path.join(params['plot_loc'],
                                         suffix+'_'+params['target']+'_signal.pdf'),
                            vmin=0,vmax=1,
                            )               
                         
    return True
# ----------------------------------------------------------------------------------------

def run_fit_continuum(params, suffix=None, start_row = None, end_row = None, 
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
	
	# First, load the datacube to be fitted.
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
    # I also need to load the SNR cube for the spaxel selection
    hdu = pyfits.open(os.path.join(params['prod_loc'],
                      '00_'+params['target']+'_snr.fits'))
    snr_cont = hdu[1].data
    hdu.close()

    # Get some info about the cube
    nrows = header1['NAXIS1']
    if not(start_row):
        start_row = 0
    if not(end_row):
	    end_row = nrows-1
	
	# Ok, what do I want to do ?
    if method == 'lowess':
        if params['verbose']:
            print '-> Starting the continuum fitting using the LOWESS approach.'
            
        fit_func = partial(brian_cof.lowess_fit, lams=lams, frac=0.05, it=5)
        # Note here the clever use of the partial function, that turns the lowess_fit
        # function from something that takes 4 arguments into something that only takes 1
        # argument ... thus perfect for the upcoming "map" functions !
        
    elif method == 'ppxf':
	    if params['verbose']:
        	print '-> Starting the continuum fitting using PPXF.' 
            #fit_func = partial(brian_ppxf.XXX)  
	        # TODO: connect to PPXF    

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
        
		# Set up the multiprocessing pool of workers
        if params['multiprocessing']:
            # Did the user specify a number of processes to use ?
            if type(params['multiprocessing']) == np.int:
                nproc = params['multiprocessing']
                pool = multiprocessing.Pool(processes = nproc, 
                                            initializer = tools.init_worker())
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(processes=None, 
                                            initializer = tools.init_worker())
			
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
                sys.stdout.flush()
			
            # Launch the fitting ! Make sure to deal with KeyBoard Interrupt properly
			# Only a problem for multiprocessing. For the rest of the code, whatever.
            try:   
                els = pool.map(fit_func, specs)
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
   	     
    print ' done !'
    
    # Ok, now, get all the fits together in one nice big cube
    # TODO: move this in it sown step ?
    cont_cube = data * 0. # That way, I conserve the possible nan's in there.

    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cube [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_'+method+'_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            conts = pickle.load(myfile)  
            myfile.close() 
        
            # Mind the shape
            cont_cube[:,:,row] = np.array(conts).T 
                   
    # Very well, now let's create a fits file to save this as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(cont_cube)
    # Make sure the WCS coordinates are included as well
    hdu1 = tools.hdu_add_wcs(hdu1,header1)
    hdu1 = tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of BRIAN is being used
    hdu1.header['BRIAN_V'] = (0.1,'brian version that created this file.')
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_'+method+'.fits')
    hdu.writeto(fn_out, clobber=True)
    
    print ' '
    
    return True
# ----------------------------------------------------------------------------------------

def run_fit_elines(params, suffix=None, start_row = None, end_row = None, 
                   ):
    ''' 
    This function fits the emission lines in the datacube, after subtracting the continuum
    derived using LOWESSS or PPXF. It is designed to use multiprocessing to speed
    things up on good computers. It deals with the data columns-per-columns, and 
    can be restarted mid-course, in case of a crash. 
    '''
    
    # Rather than launch it all at once, let's be smart in case of problems. I'll run
	# the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
	# between each row.
	
	# First, load the datacube to be fitted.
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
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
                                   '00_'+params['target']+'_snr.fits'))
    snr_cont = hdu[1].data
    hdu.close()
    
    # I also need the continuum cubes
    fn = os.path.join(params['prod_loc'],'01_'+params['target']+'_lowess.fits')
    if os.path.isfile(fn):
        hdu = pyfits.open(fn)
        cont_lowess = hdu[1].data
        hdu.close()
    fn = os.path.join(params['prod_loc'],'01_'+params['target']+'_ppxf.fits')
    if os.path.isfile(fn):
        hdu = pyfits.open(fn)
        cont_ppxf = hdu[1].data
        hdu.close()
    
    nlines = len(params['elines'].keys())
    if params['verbose']:
        print '-> Starting the emission line fitting for %2.i line(s).' % nlines 
    
    # Very well, now, perform the continuum subtraction.
    for key in params['which_cont_sub'].keys():
        if params['verbose']:
            print '   Subtracting the %s continuum from the data ...' % \
                   params['which_cont_sub'][key]
            sys.stdout.flush()
            
        llim = np.int(key.split('->')[0])
        if key.split('->')[1] == 'max':
            ulim = np.nanmax(snr_cont)+1.
        else:
        	ulim = np.int(key.split('->')[1])
         
        # Assume the continuum subtraction is error free (either lowess or models).
        if params['which_cont_sub'][key] == 'lowess':
            data[:,(snr_cont>=llim)*(snr_cont<ulim)] -= \
                 cont_lowess[:,(snr_cont>=llim) * (snr_cont<ulim)] 
        elif params['which_cont_sub'][key] == 'lowess':
            data[:,(snr_cont>=llim)*(snr_cont<ulim)] -= \
                 cont_ppfx[:,(snr_cont>=llim) * (snr_cont<ulim)]                                                      
    
    # Get some info about the cube
    nrows = header1['NAXIS1']
    if not(start_row):
        start_row = 0
    if not(end_row):
	    end_row = nrows-1 
        
    fit_func = partial(brian_elf.els_mpfit, lams=lams, be=be, params=params)
	# Note here the clever use of the partial function, that turns the els_mpfit
	# function from something that takes 4 arguments into something that only takes 1
	# argument ... thus perfect for the upcoming "map" functions !
    
    # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
	# loose the current row. Just live with it.
    for row in np.linspace(start_row,end_row, end_row-start_row+1):   
	    
	    #TODO: fit only a subset of spaxels with Halpha detected ?
		# Alright, now deal with the spaxels outside the user-chosen SNR range.
		# Replace them with nan's
        #good_spaxels = np.ones((header1['NAXIS2']))
        # if params[method+'_snr_min']:
        #     good_spaxels[snr_cont[:,row]<params[method+'_snr_min']] = np.nan
        # if params[method+'_snr_max']:
        #     good_spaxels[snr_cont[:,row]>params[method+'_snr_max']] = np.nan
		
		# Build a list of spectra to be fitted
        #specs = [data[:,i,row] * good_spaxels[i] for i in range(header1['NAXIS2'])]
        specerrs = [[data[:,i,row],error[:,i,row]] for i in range(header1['NAXIS2'])]
        
		# Set up the multiprocessing pool of workers
        if params['multiprocessing']:
            # Did the user specify a number of processes to use ?
            if type(params['multiprocessing']) == np.int:
                nproc = params['multiprocessing']
                pool = multiprocessing.Pool(processes = nproc, 
                                            initializer = tools.init_worker())
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(processes=None, 
                                            initializer = tools.init_worker())
			
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
                          suffix+'_'+params['target']+'_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')
        file = open(fn,'w')
        pickle.dump(els,file)
        file.close()
   	     
    print ' done !'
    
    return True
# ----------------------------------------------------------------------------------------

def run_make_elines_cube(params, suffix=None, prev_suffix=None):   
    ''' 
    This function is designed to construct a "usable and decent" datacube out of the
    mess generated by the emission line fitting function, i.e. out of the many pickle 
    files generated.
    '''
    
    # First, load the original datacube. I need to know how much stuff was fitted.
    hdu = pyfits.open(os.path.join(params['data_loc'],params['data_fn']))
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
    elines_params_cube = np.zeros((6*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    # And the associated errors !
    elines_params_err = np.zeros((6*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    
    if params['verbose']:
        print '-> Constructing the datacubes for the emission line fitting parameters.'
       
    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cubes [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          prev_suffix+'_'+params['target']+'_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            ms = pickle.load(myfile)  
            myfile.close() 
        
            # Get all the parameters in a big array
            ps = [item.params for item in ms]
            errs = [item.perror for item in ms]
            
            # Here, I need to make sure the errs array has a decent shape, even when the 
            # fit failed. Also store the variance = (STANDARD DEVIATION THAT COMES OUT OF 
            # MPFIT)**2 !!
            errs = [np.zeros_like(ps[0])*np.nan if item is None else item**2 for item in errs]
            
            # Fill the corresponding datacube
            elines_params_cube[:,:,row] = np.array(ps).T
            elines_params_err[:,:,row] = np.array(errs).T
            
            # now, reconstruct the full emission line spectrum
            elines_specs = np.array([brian_elf.els_spec(lams,p,be=be, 
                                                        method=params['line_profile'],
                                                        inst=params['inst']) for p in ps])
       
            # Fill the corresponding cube 
            elines_fullspec_cube[:,:,row] = elines_specs.T
            
    
    # Now, for each line, the first of these lines is the reference wavelength. 
    # Replace this by the total flux instead ! And make some plots if requested.
    for (k,key) in enumerate(params['elines'].keys()):
        
        # Calculate the sigma of the fit, in A (incl. instrument dispersion,etc ...)
        # as well as the associated error.
        zlams = elines_params_cube[6*k] * (1.+ elines_params_cube[6*k+2] / c )
        zlams_err = elines_params_cube[6*k]**2 * elines_params_err[6*k+2]/c**2
        sigma_obs_A = brian_elf.obs_sigma(elines_params_cube[6*k+3],zlams, 
                                          inst=params['inst'], 
                                          in_errs=[elines_params_err[6*k+2],zlams_err,0.] )
        
        # Compute the line flux
        elines_params_cube[6*k,:,:] = np.sqrt(2*np.pi) * elines_params_cube[6*k+1,:,:] * \
                                      sigma_obs_A[0]                           
        # What about the error ?
        elines_params_err[6*k,:,:] = 2*np.pi * (elines_params_err[6*k+1,:,:]**2 * sigma_obs_A[0] + \
                                                sigma_obs_A[1] * elines_params_cube[6*k+1,:,:]**2)           
    
    # Export the cube for each emission line parameters as a multi-extension fits file        
    # Do the same for the errors - i.e. params and errors are in two distinct cubes
    for (e,epc) in enumerate([elines_params_cube,elines_params_err]):
        hdu0 = pyfits.PrimaryHDU(None,header0)
    
        hdus = [hdu0]
        # Use the sorted keys, to ensure the same order as the fit parameters
        for (k,key) in enumerate(np.sort(params['elines'].keys())):
            hduk = pyfits.ImageHDU(epc[6*k:6*k+6,:,:])
            # Make sure the WCS coordinates are included as well
            hduk = tools.hdu_add_wcs(hduk,header1)
            # Also include a brief mention about which version of BRIAN is being used
            hduk.header['BRIAN_V'] = (0.1,'brian version that created this file.')
            # Add the line reference wavelength for future references
            hduk.header['BRIAN_L'] = (params['elines'][key][0][0], 'reference wavelength')
            hduk.header['BRIAN_C'] = ('F,I,v,sigma','Content of the cube planes')
            
            hdus.append(hduk)
            
        hdu = pyfits.HDUList(hdus=hdus)
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_elines_'+
                              ['params','perror'][e]+'.fits')
        hdu.writeto(fn_out, clobber=True)
                   
    # Very well, now let's also create a fits file to save the full emission line spectrum
    # as required.
    hdu0 = pyfits.PrimaryHDU(None,header0)
    hdu1 = pyfits.ImageHDU(elines_fullspec_cube)
    # Make sure the WCS coordinates are included as well
    hdu1 = tools.hdu_add_wcs(hdu1,header1)
    hdu1 = tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of BRIAN is being used
    hdu1.header['BRIAN_V'] = (0.1,'brian version that created this file.')
    hdu = pyfits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_elines_fullspec.fits')
    hdu.writeto(fn_out, clobber=True)
    
    print ' ' 
# ----------------------------------------------------------------------------------------

def run_plot_elines_cube(params, suffix=None, prev_suffix=None, vrange=None, 
                         sigrange=None):   
    ''' 
    This function is designed to create some plots for the emission lines, namely
    Flux, Intensity, velocity and sigma.
    '''    
    if params['verbose']:
        print '-> Making some nifty plots from the emission line fitting output.'
       
    fn = os.path.join(params['prod_loc'], prev_suffix+'_'+params['target']+
                                          '_elines_params.fits')
    
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
        for (t,typ) in enumerate(['F','I','v','sigma',]): #'h3', 'h4']):
            # Create a dedicated HDU
            tmphdu = pyfits.PrimaryHDU(this_data[t])
            # Add the WCS information
            tmphdu = tools.hdu_add_wcs(tmphdu,this_header)
            tmphdu.writeto(fn_tmp, clobber=True)
            
            # Now plot it 
            this_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                                                        '_eline-'+key+'-'+
                                                        np.str(this_header['BRIAN_L'])+
                                                        '_'+
                                                        typ+'.pdf') 
            
            # For the velocity fields, set the vmin and vmax
            if t == 0:
                my_vmin = None
                my_vmax = None
                my_cmap = None
                my_stretch = 'arcsinh'
                my_label = r'F [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [125,250,500,1000,2000,4000]
            elif t == 1:
                my_vmin = None
                my_vmax = None
                my_cmap = None
                my_stretch = 'arcsinh'
                my_label = r'I [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [50,100,200,400,800,1600]
            elif t == 2 and vrange:
                my_vmin = vrange[0]
                my_vmax = vrange[1]
                my_cmap = 'magma'
                my_stretch = 'linear'
                my_label = r'$v$ [km s$^{-1}$]'
                my_cbticks = None
            elif t ==3 and sigrange:
                my_vmin = sigrange[0]
                my_vmax = sigrange[1]
                my_cmap = 'magma'
                my_stretch = 'linear'
                my_label = r'$\sigma_{tot}$ [km s$^{-1}$]'
                my_cbticks = None
            else:
                my_vmin = None
                my_vmax = None
                my_cmap = None
                my_stretch = 'arcsinh'
                my_label = ''
                my_cbticks = None
                                                            
            brian_plots.make_2Dplot(fn_tmp,ext=0, ofn=this_ofn, contours=False, 
                                    vmin=my_vmin, vmax = my_vmax, cmap=my_cmap,
                                    stretch = my_stretch, cblabel=my_label, 
                                    cbticks = my_cbticks)                                   
    
    # Delete the temporary fits file
    os.remove(fn_tmp)
    
    # Don't forget to close the initial hdu ...
    hdu.close()
    
    return True
# ----------------------------------------------------------------------------------------

def run_find_structures(params, suffix=None, interactive_mode=True):   
    ''' 
    This function is designed to identify structures (e.g. HII regions) in the data from
    a 2D image (i.e. an line intensity map), and save them to a pickle file. When
    interactive_mode=True, the user can manually refine the selection.
    '''
    
    if params['verbose']:
        print '-> Starting the semi-automated procedure for structure identification.'
    
    
    # Where am I going to save the apertures information ?
    fn_ap = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_aplist.pkl')
    
    # Do we want to build apertures based on multiple maps ? Loop, one after the other.
    for key in params['ap_map_lines']:
        print '   Loading eline %s' %key 
        # Is the file already there ? Do I want to overwrite it ? And/or use its content ?
        if os.path.isfile(fn_ap):
            print '    '
            print '   Existing aperture list (%s)' % fn_ap.split('/')[-1]
            print '   Type [a] to load this aperture list and edit it manually'
            print '        [b] to start from scratch, detecting structures automatically'
            print '        [c] to do nothing, and continue to the next processing step'
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
            return True
    
        # Now, open the elines param datacube, and extract the Flux map I want to detect
        # stuctures from.
        fn = os.path.join(params['prod_loc'],'04_'+params['target']+'_elines_params.fits')
        hdu = pyfits.open(fn)
        header0 = hdu[0].header
        plane = np.sort(params['elines'].keys()).tolist().index(key)
        data = hdu[plane+1].data[0]
        hdu.close()
        
        # Launch the aperture finding routine
        apertures = brian_plots.build_ap_list(data, start_aps = start_aps, 
                                              radius = params['ap_radius'])    
        
        # Only if the user wants to save the apertures, do it
        if apertures:
            # Save the results for later use
            f = open(fn_ap, 'w')
            pickle.dump(apertures,f)
            f.close()
        
    return True
        
             
    