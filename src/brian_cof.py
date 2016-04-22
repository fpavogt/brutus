# -*- coding: utf-8 -*-
#
# This file contains functions related to the continuum fitting inside BRIAN, except for
# anything related to PPXF. 
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

import brian_tools as tools

# ----------------------------------------------------------------------------------------

def lowess_fit(spec, lams, frac=0.05, it=5):
    '''
    This function fits a spectrum using a LOWESS (Locally Weighted Scatterplot Smoothing)
    technique, described in: 
    Cleveland, W.S. (1979) Robust Locally Weighted Regression and Smoothing Scatterplots. 
    Journal of the American Statistical Association 74 (368): 829-836.
    
    This is robust to outliers (hot pixels, cosmics), and is also efficient to ignore 
    emission lines. frac=0.05 and it=5 seem to work very fine for spectra of any SNR, both
    lousy with no continuum, and good ones in the center of galaxies - modulo the stellar
    absorption features which are of course "ignored" by the LOWESS routine.
    '''
    
    fit = lowess(spec,lams,frac=0.05, it=5, is_sorted=True, missing = 'drop', 
                 return_sorted=False)
	                                               
    return fit
# ----------------------------------------------------------------------------------------      

	