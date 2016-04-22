# -*- coding: utf-8 -*-
#
# This file contains tools for the BRIAN routines to fit the stellar continuum and the 
# emission lines in an IFU data cube (i.e. MUSE).
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
import signal 
 
 # ----------------------------------------------------------------------------------------      
  
def init_worker():
    '''
    THis function is designed to handle KeyboardInterrupt during multiprocessing.
    See https://noswap.com/blog/python-multiprocessing-keyboardinterrupt
    '''
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    
# ----------------------------------------------------------------------------------------      
  
def hdu_add_wcs(newheader,refheader):
    '''
    This function adds the WCS coordinates from a reference header to a new header.
    '''
    
    newheader.header['CRPIX1'] = refheader['CRPIX1']
    newheader.header['CD1_1'] = refheader['CD1_1']
    newheader.header['CTYPE1'] = refheader['CTYPE1']
    newheader.header['CUNIT1'] = refheader['CUNIT1']
    newheader.header['CRPIX2'] = refheader['CRPIX2']
    newheader.header['CD2_2'] = refheader['CD2_2']
    newheader.header['CTYPE2'] = refheader['CTYPE2']
    newheader.header['CUNIT2'] = refheader['CUNIT2']
    newheader.header['CD1_2'] = refheader['CD1_2']
    newheader.header['CD2_1'] = refheader['CD2_1']
    newheader.header['CRVAL1'] = refheader['CRVAL1']
    newheader.header['CRVAL2'] = refheader['CRVAL2']
    
    return newheader
# ----------------------------------------------------------------------------------------      
    
def hdu_add_lams(newheader,refheader):
    '''
    This function adds the Wavelength information from a reference header to a new header.
    '''
    
    newheader.header['CTYPE3'] = refheader['CTYPE3'] 
    newheader.header['CUNIT3'] = refheader['CUNIT3']   
    newheader.header['CD3_3'] = refheader['CD3_3']
    newheader.header['CRPIX3'] = refheader['CRPIX3']
    newheader.header['CRVAL3'] = refheader['CRVAL3']
    newheader.header['CD1_3'] = refheader['CD1_3']
    newheader.header['CD2_3'] = refheader['CD2_3']
    newheader.header['CD3_1'] = refheader['CD3_1']
    newheader.header['CD3_2'] = refheader['CD3_2']
 
    return newheader    