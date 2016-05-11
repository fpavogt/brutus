# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains several function and tools used by the brutus routines to correct
 emission line fluxes for galactic and extragalactic reddening.

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import scipy as sp
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from brutus_metadata import *

# ----------------------------------------------------------------------------------------      

def D_fm90(x,gamma,x0):
    ''' The UV points from Fitzpatrick & Massa (1990). 
    
    :Args:
        x: float
           Wavelength in 1/microns.
        gamma: float
               The gamma factor.
        x0: float
            The x0 factor.
    '''
    return x**2/((x**2-x0**2)**2+x**2*gamma**2)
# ---------------------------------------------------------------------------------------- 
     
def F_fm90(x):
    '''The UV points from Fitzpatrick & Massa (1990). 
        
    :Args:
        x: float
           Wavelength in 1/microns.
    
    '''
    out = np.zeros_like(x)
    for (i,item) in enumerate(x):
        if item >= 5.9:    
            out[i] = 0.5392*(item-5.9)**2 + 0.05644*(item-5.9)**3
        else :
            out[i] = 0.0
    return out
# ----------------------------------------------------------------------------------------  
            
def f99_alebv(lams,rv=3.1): # lams in Angstroem
    '''The spline-interpolated Fitzpatrick, PASP (1999) Alambda/E(B-V) function.
    
    :Args: 
        lams: 1-D numpy array of floats
              The wavelengths in Angstroem.
        rv: float [default: 3.1]
            The Rv value.
    
    '''
    
    my_lams = 1.0e4/lams 
     
    # Make sure this is an array, or the code will crash
    try :
        len(my_lams)
    except :
        my_lams = np.array([my_lams])
       
    # First, the anchor points (in 1/lambda)
    uv_anchors = np.array([3.704,3.846])
    optical_anchors = np.array([1.667, 1.828, 2.141, 2.433])
    ir_anchors = np.array([0.0,0.377, 0.820])    
             
    # ---- The parameters for the UV points ----
    c2 = -0.824 + 4.717/rv
    c1 = 2.030 - 3.007 * c2
    x0 = 4.596
    gamma = 0.99
    c3 = 3.23
    c4 = 0.41
    D = D_fm90(uv_anchors,gamma,x0)
    F = F_fm90(uv_anchors)    
    kx_V = c1 + c2 * uv_anchors + c3 * D + c4 * F
    uv_spline = kx_V + rv
    
    # ---- For the optical anchors ----
    optical_spline = np.zeros_like(optical_anchors)
    optical_spline[0] = -0.426 + 1.0044 * rv
    optical_spline[1] = -0.050 + 1.0016 * rv
    optical_spline[2] =  0.701 + 1.0016 * rv 
    optical_spline[3] = +1.208 + 1.0032 * rv - 0.00033 * rv**2 
    # WARNING : typo in F99 ?
    # See their Table 4 ... with '-1.208', it doesn't work !

    # ---- For the IR anchors ----
    ir_spline = rv/3.1 * np.array([0.0,0.265,0.829])
    
    # Ok, now, perform the spline interpolation
    wave = np.append(np.append(ir_anchors, optical_anchors), uv_anchors)
    alebv = np.append(np.append(ir_spline, optical_spline), uv_spline)   
    interp_func = interpolate.interp1d(wave, alebv, kind='cubic', 
                                       bounds_error = False) 
    # bounds_error = False -> Will return a NaN is outside the allowed range 
    
    # Finally, compute the value at the wavelengths of interest
    out = np.zeros_like(my_lams)
    
    # If lam > 2700, use the spline function
    out[lams>2700] = interp_func(my_lams[lams>2700])
    # If lam <= 2700, use the fm90 function
    Dl = D_fm90(my_lams[lams<=2700], gamma, x0)
    Fl = F_fm90(my_lams[lams<=2700])
    out[lams<=2700] = rv + c1 + c2 * my_lams[lams<=2700] + c3 * Dl + c4 * Fl
    
    return out  
# ----------------------------------------------------------------------------------------  
# Get E_lambda-V/E_B-V following Fischera & Dopita05. Unlike Vogt et al. (2013), follow 
# the paper more closely, and allow for varying values of Rv using polynomials given in 
# the text. In other words, use BOTH rv and rva to characterize the attenuation curve.
def fd05_elvebv(lams, rv = 3.08, rva = 4.3): # lams in Angstroem
    '''Compute E(lambda-V)/E(B-V) according to Fischera&Dopita (2005), given Rv and Rva.
    
    This function reproduces Table 1 in the article. Unlike in Appendix 1 of Vogt et al.
    (2013), this function relies on the polynomials defined in the text to reproduce the 
    values of E(lambda-V)/E(B-V) for any Rv and Rva.
    
    :Args: 
        lams: 1-D numpy array of float
              The wavelengths at which to evaluate the attenutation.
        rv: float [default: 3.08]
            The Rv value, following Fischera & Dopita.
        rva: float [default: 4.5]
             The Rva value, following Fischera & Dopita.
    
    :Returns:
        out: 1-D numpy aarray of float
             The value of E(lambda-V)/E(B-V), as in Table 1 and Fig.3 of the article.
             
    :Notes:
        By default, Rv = 3.08 and Rva = 4.3, which gives a curve vey close to that of 
        Calzetti (2000) below the 2000 Ansgtroem bump. Thanks to R. Sutherland at ANU/RSAA
        for sharing his reddening code (included in MAPPINGS V), which helped me implement 
        this function properly.
    '''
    
    # Make sure this is an array (for the overall consistency)
    try :
        len(lams)
    except :
        lams = np.array([lams])
        
    # Rv < 3.08 is probably a bad thing ...    
    if rv < 3.08:
        sys.exit(" Rv <3.08 not supported !")
    
    # The wavelength nodes from Table 1 in microns
    lam_nodes = np.array([0.0912, 0.1053, 0.1216, 0.1550, 0.1906, 0.2175, 0.2480, 0.3000,
                          0.3650, 0.4400, 0.5480, 0.7200, 1.0300, 1.2390, 1.6490, 2.1920,
                          3.5920, 4.7770])
    
    # Do it all for the Av=1 curves. According to the paper:
    #  If the approximation obtained for AV= 1 is used for a turbulent medium where 
    # AV= 10 the accuracy is better than ∼5%.
    
    # The curve for Rv=Rva=3.08
    elv_ebv308a1 = np.array([12.066, 9.694, 7.328, 5.018, 5.214, 6.937, 4.573, 2.967, 
                             1.924, 1.000, 0.000,-0.991,-1.871,-2.171,-2.498,-2.716,
                             -2.932, -2.992])
                              
    # The x from Eqt. 8. Al/Av = (Elv/Ebv*1/Rv) + 1
    xs = elv_ebv308a1/3.08 + 1.0
    
    # The a coefficients, as polynomial using the b coefficients from Table 2.
    fa1 = np.poly1d([0.0, 1.64220, -3.62163, 3.42939, 0.00495676])
    fa2 = np.poly1d([-2.199, 5.11725, -0.668874, -0.708645, 0.0393932])
    fa3 = np.poly1d([21.5465, -19.2137, 5.17489, -0.484444, 0.0174931])
    
    # By definition, following Par.2 in Sec 3.1.:
    alphav = 3.08 
    
    # Eqt. 11
    m = 1 + 0.275 * (rv - alphav)
    
    # Eqt. 10
    alphava = m**-1*(rva - rv) + alphav
    
    # Eqt. 9
    z = (alphava-1.0)**-1
    a1 = fa1(z)
    a2 = fa2(z)
    a3 = fa3(z)
    
    # Eqt. 8
    alav = a1 * np.log10(xs)
    alav += a2 * (np.log10(xs))**2
    alav += a3 * (np.log10(xs))**3
    alav = 10.**alav
    
    # Reproduce the values of table 1. These are approximated using the polynomials, and 
    # according to the article (for the table, set rv=3.08 and rva =...):
    # If applied to situations where the absolute attenuation AV is either 1 or 10, the 
    # approximations are better than 1.3%.  In general the accuracy is somewhat lower.  
    # If the approximation obtained for AV= 1 is used for a turbulent medium where AV= 10 
    # the accuracy is better than ∼5%. The largest errors occurs at low optical depths. 
    # For x≥1 the approximation is better than ∼2%.
    elbebv = (alav-1.) * rva
    
    # Very well, now, I need to interpolate between these nodes, to get a smooth function 
    # at "any" wavelength. Use Akima splines. Fit the curve in 1/lam space, because the
    # points are a lot smoother then. Reverse their order, because x must be strictly 
    # ascending. To be clear, this fits elb_ebv as a function of 1/lambda in 1/microns
    akima_spline = sp.interpolate.Akima1DInterpolator(1./lam_nodes[::-1],elbebv[::-1]) 
    
    return akima_spline(1./(lams*1e-4))
# ----------------------------------------------------------------------------------------  

def cal00_ke(lams,rv):
    '''The Calzetti (2000) law. lams in Angstroem. 
    
    :Notes:
        This does NOT include the 0.44 correction factor to apply for stellar extinction.
    '''
    
    # Make sure lams is an array (for the overall consistency)
    try :
        len(lams)
    except :
        lams = np.array([lams])
    
    my_lams = 1.0e4/lams             
    ke = np.zeros_like(lams,dtype='d')
    for i in range(len(lams)):
        if lams[i]/1.0e4 < 0.63 and lams[i]/1.0e4 >= 0.12 :
            ke[i] = 2.659 * (-2.156 +1.509 * my_lams[i] - 0.198*my_lams[i]**2 +
                             0.011 * my_lams[i]**3) + rv
        elif lams[i]/1.0e4 < 2.20 and lams[i]/1.0e4 >= 0.63 :
            ke[i] = 2.659 * (-1.857 +1.040 *my_lams[i]) + rv
    return ke 
# ---------------------------------------------------------------------------------------- 

# Following Cardelli, Clayton & Mathis (1989)
def ccm89_alav(lams, rv): # lams in Angstroem
    '''
    The law from Caredelli, Clayton and Mathis (1989). lams in Angstroem.
    
    '''
    
    # Make sure lams is an array 
    try :
        len(lams)
    except :
        lams = np.array([lams])

    my_lams = 1.0e4/lams
    
    a = np.zeros_like(lams, dtype ='d')
    b = np.zeros_like(lams, dtype ='d')
    alav = np.zeros_like(lams, dtype ='d')
    
    for i in range(len(my_lams)):
        x = my_lams[i]
        if 0.3 <= x and x < 1.1 :
            a[i] =  0.574 * x**1.61
            b[i] = -0.527 * x**1.61
        elif 1.1 <= x and x < 3.3 :
            y = x-1.82
            a[i] = 1. + 0.17699 * y -0.50447 * y**2 - 0.02427 * y**3 +\
                0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
            b[i] = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 -\
                0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
        elif 3.3 <= x and x<= 8 :
            if x < 5.9 : 
                fa = 0
                fb = 0
            else :
                fa = -0.04473 * (x-5.9)**2 - 0.009779 * (x-5.9)**3
                fb =  0.2130 * (x-5.9)**2 + 0.1207 * (x-5.9)**3
            a[i] =  1.752 - 0.316 * x - 0.104/((x-4.67)**2 + 0.341) + fa
            b[i] = -3.090 + 1.825 * x + 1.206/((x-4.62)**2 + 0.263) + fb
    
    alav = a + b/rv

    return alav
# ---------------------------------------------------------------------------------------- 

def alam(lams, ab, av, curve='f99', rv=None, rva=4.5):
    '''Calculate the value of A_lambda for different extinction/attenuation laws.
    
    :Args:
        lams: 1-D numpy array of float
              The values at which to evaluate the attenuation/extinction.
        ab: float
            The extinction in B (e.g. according to NED)
        av: float
            The extinction in V (e.g. according to NED) 
        curve: string [default: 'f99']
               Which extinction/attenuation curve to use.
        rv: float
            The value of Rv.
        rva: float [default: 4.5]
            The value of Rva, if using curve='fd05'.
        
    :Notes: 
        To reproduce the values given in NED, set curve='f99', and rv = 3.1. The default 
        value of Rv will depend on the curve used (if not specified). rv=3.1 for 'f99', 
        rv=3.08 for 'fd05', rv=4.05 for 'cal00' and 'cal00*', and rv=3.1 for 'ccm89'. 
        'cal00*' includes a 0.44 correction factor to derive the stellar extinction.  
        
    '''
    if curve == 'f99':
        if rv is None:
            rv = 3.1
        out = f99_alebv(lams,rv) * (ab-av) 
        return out
        
    if curve == 'fd05':
        # From Fischera& Dopita (2005), see also Vogt et al. (2013).
        if rv is None:
            rv = 3.08   
        red = fd05_elvebv(lams, rv=rv, rva=rva)
        return (red + rva) * (ab-av)
        
    elif 'cal00' in curve:
        if rv is None :
            rv = 4.05
        corr = 1.0
        if curve == 'cal00*':
            corr = 0.44
        return cal00_ke(lams,rv)*(ab-av)*corr
        
    elif curve =='ccm89':
        if rv is None:
            rv = 3.1
        return ccm89_alav(lams,rv)*av      
# ---------------------------------------------------------------------------------------- 

def galactic_red(lams, ab, av, curve='f99', rv = 3.1, rva = 4.3):
    '''Calculate the galactic extinction for a given sight-line, using Ab and Av.
    
    :Args:
        lams: 1-D numpy array of float
             The values at which to evaluate the attenuation/extinction.
        ab: float
            The extinction in B (e.g. according to NED).
        av: float
            The extinction in V (e.g. according to NED) .
        curve: string [default: 'f99']
               Which extinction/attenuation curve to use.
        rv: float
            Rv value.
        rva: float [default: 4.5]
             The value of Rva, if using curve='fd05'.
         
    :Returns:
        out: 1-dNumpy array
             The ratio of Flux_(unreddened)/ Flux_(observed).
    :Notes:
        To follow NED, set curve='f99', and rv = 3.1. The default 
        value of Rv will depend on the curve used (if not specified). rv=3.1 for 'f99', 
        rv=3.08 for 'fd05', rv=4.05 for 'cal00' and 'cal00*', and rv=3.1 for 'ccm89'. 
        'cal00*' includes a 0.44 correction factor to derive the stellar extinction.
        
    '''
    tau = alam(lams,ab,av,curve=curve,rv=rv) * (2.5*np.log10(np.e))**-1
          
    return np.exp(tau)
# ---------------------------------------------------------------------------------------- 

def extragalactic_red(lam, hahb, hahb_0, curve = 'fd05', rv = None, rva=4.3):
    '''Calculate the extragalactic attenuation (based on Halpha/Hbeta).
    
    :Args:
        lam: float
              The wavelength at which to evaluate the attenuation/extinction.
        hahb: float or numpy array
              The observed Halpha/Hbeta flux ratio. Can be a 2-D map.
        hahb_0: float
                The reference theoretical unreddened Halpha/Hbeta flux ratio.
        curve: string [default: 'fd05']
               Which extinction/attenuation curve to use.
        rv: float [default: None]
            Rv value.
        rva: float [default: 4.3]
             The value of Rva, if using curve='fd05'. 
    :Returns:
        out: 1-D numpy array
             The ratio of Flux_(unattenuated)/ Flux_(attenuated)
    :Notes:
        See also Vogt+ (2013), Appendix A. By default, Rv = 3.08 and Rva = 4.3, which 
        gives a curve vey close to that of Calzetti (2000) below the 2000 Ansgtroem bump. 
        Supported curves are 'fd05', 'f99','ccm89','cal00' and 'cal00*'.'cal00*' includes 
        a 0.44 correction factor to derive the stellar extinction.
        
    '''
    if curve == 'fd05' :
        elvebv = fd05_elvebv(lam)
        ehavebv = fd05_elvebv(ha)	
        ehbvebv = fd05_elvebv(hb) 
        out = (hahb/hahb_0)**(-(elvebv + rva)/(ehavebv-ehbvebv))  
        return out
    if curve == 'ccm89':
        if rv is None:
            rv = 3.1
        alav = ccm89_alav(lam, rv)
        aHaav = ccm89_alav(ha, rv)
        aHbav = ccm89_alav(hb, rv)
        return (hahb/hahb_0)**(-alav/(aHaav-aHbav))
        
    elif 'cal00' in curve:
        if rv == None:
            rv=4.05
        ke = cal00_ke(lam,rv)
        keHa = cal00_ke(ha,rv)
        keHb = cal00_ke(hb,rv)
        if curve == 'cal00*':
            ke *=0.44
            keHa *= 0.44
            keHb *= 0.44
        return (hahb/hahb_0)**(-ke/(keHa-keHb))
        
    elif curve == 'f99':
        if rv == None:
            rv=3.1
        alebv = f99_alebv(lam,rv)
        aHaebv = f99_alebv(ha,rv)
        aHbebv = f99_alebv(hb,rv)
        return (hahb/hahb_0)**(-alebv/(aHaebv-aHbebv))
# ----------------------------------------------------------------------------------------         

def hahb_to_av(hahb,hahb_0, curve = 'fd05', rv = None, rva=4.3):
    '''A function to get the reddening in Av mag, from Ha/Hb ratio.
    
    :Args: 
        ha_hb: 1-D numpy array
               The flux ratio of Halpha to Hbeta (NOT the log!)
        ha_hb_0: float
                 The theoretical value of Halpha/Hbeta. (2.86 for SB, more for AGNs)
        curve: string [default: 'fd05']
               Which reddeing law to use ?
        rv: float [default: None]
            The value of Rv to use. None for the default.
    
    :Returns: 
        out: 1-D numpy array
             The corrseponding values of Av.
    :Notes:
        See also Vogt+ (2013), Appendix A. By default, Rv = 3.08 and Rva = 4.3, which 
        gives a curve vey close to that of Calzetti (2000) below the 2000 Ansgtroem bump.
        
    '''
    if curve == 'fd05' :
        ehavebv = fd05_elvebv(ha)	
        ehbvebv = fd05_elvebv(hb) 
        out = -2.5*np.log10(hahb/hahb_0)*(rva/(ehavebv-ehbvebv))  
        return out
    if curve == 'ccm89':
        if rv == None:
            rv = 3.1
        aHaav = ccm89_alav(ha, rv)
        aHbav = ccm89_alav(hb, rv)
        return -2.5*np.log10(hahb/hahb_0)*(1./(aHaav-aHbav))
        
    elif curve == 'cal00':
        if rv == None:
            rv=4.05
        keHa = cal00_ke(ha,rv)
        keHb = cal00_ke(hb,rv)
        return -2.5*np.log10(hahb/hahb_0)*(rv/(keHa-keHb))
        
    elif curve == 'f99':
        if rv == None:
            rv=3.1
        aHaebv = f99_alebv(ha,rv)
        aHbebv = f99_alebv(hb,rv)
        return -2.5*np.log10(hahb/hahb_0)*(rv/(aHaebv-aHbebv))
# ---------------------------------------------------------------------------------------- 

def check():
    '''
    Make some quick plots to test the reddening functions of brutus, by comparing with
    the original papers. And make sure I did not mess things up ...
    '''
    
    # 1) Check Cardelli89 and Fitzpatrick99: reproduce Fig. 7 of Fitzpatrick99    
    x = np.arange(0.0001,8.75,0.001)
    lams = 1/x*1.e4
     
    plt.close(97)
    plt.figure(97, figsize = (15,7))
    gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.1,right=0.95,bottom=0.15,top=0.9,wspace=0.05,hspace=0.05 )
     
    ax1 = plt.subplot(gs[0,0])
     
    ls = ['r','m','g','b']
     
    for (i,r) in enumerate([2.3,3.1,4.0,5.5]):
        f99_curve = f99_alebv(lams,r) - r
        ccm89_curve = (ccm89_alav(lams,r)*r) -r
        
        ax1.plot(x,f99_curve,ls[i]+'-', label=r'R$_V$ = '+np.str(r), linewidth=2)
        ax1.plot(x,ccm89_curve,'k--')

    ax1.grid(True)
    ax1.set_xlabel(r'1/$\lambda$ [$\mu$m$^{-1}$]')
    ax1.set_ylabel(r'E($\lambda$-V)/E(B-V)')
    ax1.legend(loc='lower right',fontsize=15)
    ax1.text(0.1,11,'Fitzpatrick (1999) - Fig. 7; F99 (-) vs CCM89 (- -)')
    ax1.set_xlim([0,10])
    ax1.set_ylim([-6,13])
    # Add the MUSE range because I can
    ax1.axvline(x=1./(4750./1e4),ymax=0.65,c='k',ls='-', lw=2)
    ax1.axvline(x=1./(9350./1e4),ymax=0.65,c='k',ls='-', lw=2)
    ax1.text(1.7,7,r'MUSE range (4750\AA $\rightarrow$ 9350\AA)',ha='center')
    
    plt.show()
    
    # --------------------------------------------------------------------------
    # 2) Compare with Fischera & Dopita 2005  
    
    x = np.arange(0.0001,7.5,0.1)
    lams = 1/x*1.e4
     
    plt.close(98)
    plt.figure(98, figsize=(10,8))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,0.5], width_ratios=[1])
    gs.update(left=0.15,right=0.9,bottom=0.1,top=0.95,wspace=0.05,hspace=0.05 )
    
    r = 3.1
    f99_curve = f99_alebv(lams,r) - r
    ccm89_curve = (ccm89_alav(lams,r)*r) -r
    cal00_curve = cal00_ke(lams,4.05) - 4.05   
    fd05_curve = fd05_elvebv(lams,rv=3.08,rva=4.3)
    
    ax1 = plt.subplot(gs[0,0])
    ax1.plot(x,f99_curve,'r-', linewidth=1,label=r'F99, R$_V$='+np.str(r))
    ax1.plot(x,ccm89_curve,'g-',label=r'CCM89, R$_V$='+np.str(r))
    ax1.plot(x,cal00_curve,'b-',label=r'Cal00, R$_V$=4.05')
    ax1.plot(x,fd05_curve,'k--',linewidth=2, label='FD05,R$_V$=3.08, R$_V^{A}$=4.3')
    ax1.grid(True)
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'E($\lambda$-V)/E(B-V)')
    ax1.legend(loc='lower right',fontsize=15)
    ax1.set_xlim([0.5,4.5])
    ax1.set_ylim([-5,5])
    ax1.axvline(x=1./(4750./1e4),ymax=0.65,c='k',ls='-', lw=2)
    ax1.axvline(x=1./(9350./1e4),ymax=0.65,c='k',ls='-', lw=2)
    ax1.text(1.7,2.5,r'MUSE range (4750\AA $\rightarrow$ 9350\AA)',ha='center')
    
    ax2 = plt.subplot(gs[1,0])
    ax2.axhline(y=0,color =(0.4,0.4,0.4), linewidth=2)
    ax2.plot(x,f99_curve-fd05_curve,'r-', linewidth=1,label=r'F99, R$_V$='+np.str(r))
    ax2.plot(x,ccm89_curve-fd05_curve,'g-',label=r'CCM89, R$_V$='+np.str(r))
    ax2.plot(x,cal00_curve-fd05_curve,'b-',label=r'Cal00, R$_V$='+np.str(r))    
    ax2.grid(True)
    ax2.set_ylabel(r'$\Delta$E($\lambda$-V)/E(B-V)')
    #ax2.legend(loc='lower right')
    ax2.set_xlim([0.5,4.5])
    #ax1.ylim([-5,5])
    ax2.set_xlabel(r'1/$\lambda$ [$\mu$m$^{-1}$]')
    ax2.axvline(x=1./(4750./1e4),c='k',ls='-', lw=2)
    ax2.axvline(x=1./(9350./1e4),c='k',ls='-', lw=2)
       
    plt.show()
    
    # --------------------------------------------------------------------------
    # 3) Reproduce Fig.1 of Calzetti (2001) 
    
    x = np.arange(0.0001,7.5,0.01)
    lams = 1/x*1.e4
     
    plt.close(99)
    plt.figure(99, figsize=(10,8))
    gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.15,right=0.9,bottom=0.1,top=0.95,wspace=0.05,hspace=0.05 )
    
    f99_2 = f99_alebv(lams,2.0)
    f99_31 = f99_alebv(lams,3.1)
    f99_5 = f99_alebv(lams,5.0)
    ccm89_2 = ccm89_alav(lams,2.0)*2.0
    ccm89_31 = ccm89_alav(lams,3.1)*3.1
    ccm89_5 = ccm89_alav(lams,5.0)*5.0
    #ccm89_curve = (ccm89_alav(lams,r)*r) -r
    cal00_curve = cal00_ke(lams,4.05)   
    fd05_curve = fd05_elvebv(lams,rv=3.08,rva=4.3)+4.3
    
    ax1 = plt.subplot(gs[0,0])
    ax1.plot(np.log10(lams/1.e4),f99_2,'g:', linewidth=1,label=r'F99, R$_V$=2.0')
    ax1.plot(np.log10(lams/1.e4),f99_31,'g-', linewidth=1,label=r'F99, R$_V$=3.1')
    ax1.plot(np.log10(lams/1.e4),f99_5,'g--', linewidth=1,label=r'F99, R$_V$=5.0')
    ax1.plot(np.log10(lams/1.e4),ccm89_2,'r:', linewidth=1,label=r'CCM89, R$_V$=2.0')
    ax1.plot(np.log10(lams/1.e4),ccm89_31,'r-', linewidth=1,label=r'CCM89, R$_V$=3.1')
    ax1.plot(np.log10(lams/1.e4),ccm89_5,'r--', linewidth=1,label=r'CCM89, R$_V$=5.0')
    ax1.plot(np.log10(lams/1.e4),cal00_curve,'b-',label=r'Cal00, R$_V$=4.05')
    ax1.plot(np.log10(lams/1.e4),fd05_curve,'k--',linewidth=2, label='FD05,R$_V$=3.08, R$_{V}^{A}$=4.3')

    ax1.grid(True)
    ax1.set_ylabel(r'A$_\lambda$/E(B-V)')
    ax1.set_xlabel(r'log[$\lambda$($\mu$m)]', labelpad=10)
    ax1.legend(loc='upper right', fontsize=15)
    ax1.set_xlim([-1,0.3])
    ax1.set_ylim([0,15])
    
    # Add the MUSE range because I can
    ax1.axvline(x=np.log10(4750./1e4),ymax=0.5,c='k',ls='-', lw=2)
    ax1.axvline(x=np.log10(9350./1e4),ymax=0.5,c='k',ls='-', lw=2)
    ax1.text(-0.15,8.1,r'MUSE range (4750\AA $\rightarrow$ 9350\AA)',ha='center')
    
       
    plt.show()
# ----------------------------------------------------------------------------------------


