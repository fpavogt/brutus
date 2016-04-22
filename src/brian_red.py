# -*- coding: utf-8 -*-
#
# This file contains several function and tools used by the BRIAN routines to correct
# emission line fluxes for galactic and extragalactic reddening.
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from brian_metadata import *

# ----------------------------------------------------------------------------------------      

def D_fm90(x,gamma,x0):
    ''' 
    The UV points from Fitzpatrick & Massa (1990). x is in 1/microns
    '''
    return x**2/((x**2-x0**2)**2+x**2*gamma**2)
# ---------------------------------------------------------------------------------------- 
     
def F_fm90(x):
    '''
    The UV points from Fitzpatrick & Massa (1990). x is in 1/microns
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
    '''
    The spline-interpolated Fitzpatrick, PASP (1999) Alambda/E(B-V) function.
    lams in Ansgtroem
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
# TODO
# Get E_lambda-V/E_B-V following Fischera & Dopita05, Vogt et al. (2013)
'''
def fd05_elvebv(lams): # lams in Angstroem
# Make sure this is an array (for the overall consistency)
    try :
        len(lams)
    except :
        lams = np.array([lams])
    
    my_lams = 1.0e4/lams    
    out = -4.61777 + 1.41612 * my_lams + 1.52077 * my_lams**2 \
        - 0.63269 * my_lams**3 +0.07386 * my_lams**4    
      # WARNING : this is ONLY valid over the range 3000<lambda<12000 !
      # Make it return 'NaN' elsewhere    
    out[lams>12390] = np.nan
    out[lams<2480] = np.nan
    
    return out
'''
# ----------------------------------------------------------------------------------------  

def cal00_ke(lams,rv):
    '''
    The Calzetti (2000) law. lams in Angstroem.
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
    The law from Caredelli, CLayton and Mathis (1989). lams in Angstroem.
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

def alam(lams,ab,av,curve='f99',rv=None, cal00_stellar=True):
    '''
    Calculate the value of A_lambda for different curves -> bring everyone together !
    Input : lambda(angstroem), ab & av (e.g. from NED), rv and extinction curve
            cal00_stellar (True/False), whethert o add the 0.44 correction factor for 
            stellar extinction (vs gas)
    '''
    if curve == 'f99':
        if rv == None:
            rv = 3.1
        out = f99_alebv(lams,rv) * (ab-av) 
        return out
        
    #if curve == 'fd05':
    #    # From Fischera& Dopita (2005), see also Vogt et al. (2013)
    #    rv = 4.5     
    #    red = fd05_elvebv(lams)
    #    return (red + rv)*(ab-av)
    elif curve == 'cal00':
        if rv == None :
            rv = 4.05
        corr = 1.0
        if cal00_stellar:
            corr = 0.44
        return cal00_ke(lams,rv)*(ab-av)*corr
        
    elif curve =='ccm89':
        if rv == None:
            rv = 3.1
        return ccm89_alav(lams,rv)*av      
# ---------------------------------------------------------------------------------------- 

def galactic_red(lams,ab,av,curve='f99',rv = 3.1, cal00_stellar = False):
    '''
    Calculate galactic reddening (flux_ratio). 
    Input: curve (which function) and lams, the wavelength in Angstroem,
    ab and av, the extinction from NED
    Return F unreddened/ F observed
    '''
    tau = alam(lams,ab,av,curve=curve,rv=rv, cal00_stellar=cal00_stellar) * \
          (2.5*np.log10(np.e))**-1
          
    return np.exp(tau)
# ---------------------------------------------------------------------------------------- 

def extragalactic_red(lam,ha_hb,ha_hb_0, curve = 'fd05', rv = None):
    '''
    Calculate extragalactic reddening (based on Halpha/Hbeta) 
    Input: wavelengths (angstroem), Halpha/Hbeta (observed), 
    Halpha/Hbeta (reference), extinction curve, rv
    Return F unreddened/ F observed
    '''
    #if curve == 'fd05' :
    #    elvebv = fd05_elvebv(lam)
    #    ehavebv = fd05_elvebv(halpha)	
    #    ehbvebv = fd05_elvebv(hbeta) 
    #    out = (ha_hb/ha_hb_0)**(-(elvebv+4.5)/(ehavebv-ehbvebv))  
    #    return out
    if curve == 'ccm89':
        if rv == None:
            rv = 3.1
        alav = ccm89_alav(lam, rv)
        aHaav = ccm89_alav(ha, rv)
        aHbav = ccm89_alav(hb, rv)
        return (ha_hb/ha_hb_0)**(-alav/(aHaav-aHbav))
        
    elif curve == 'cal00':
        if rv == None:
            rv=4.05
        ke = cal00_ke(lam,rv)
        keHa = cal00_ke(halpha,rv)
        keHb = cal00_ke(hbeta,rv)
        return (ha_hb/ha_hb_0)**(-ke/(keHa-keHb))
        
    elif curve == 'f99':
        if rv == None:
            rv=3.1
        alebv = f99_alebv(lam,rv)
        aHaebv = f99_alebv(halpha,rv)
        aHbebv = f99_alebv(hbeta,rv)
        return (ha_hb/ha_hb_0)**(-alebv/(aHaebv-aHbebv))
# ----------------------------------------------------------------------------------------         

def hahb_to_av(ha_hb,ha_hb_0, curve = 'fd05', rv = None):
    '''
    A function to get the reddening in Av mag, from Ha/Hb ratio
    '''
    #if curve == 'fd05' :
    #    ehavebv = fd05_elvebv(halpha)	
    #    ehbvebv = fd05_elvebv(hbeta) 
    #    out = -2.5*np.log10(ha_hb/ha_hb_0)*(4.5/(ehavebv-ehbvebv))  
    #    return out
    if curve == 'ccm89':
        if rv == None:
            rv = 3.1
        aHaav = ccm89_alav(ha, rv)
        aHbav = ccm89_alav(hb, rv)
        return -2.5*np.log10(ha_hb/ha_hb_0)*(1./(aHaav-aHbav))
        
    elif curve == 'cal00':
        if rv == None:
            rv=4.05
        keHa = cal00_ke(halpha,rv)
        keHb = cal00_ke(hbeta,rv)
        return -2.5*np.log10(ha_hb/ha_hb_0)*(rv/(keHa-keHb))
        
    elif curve == 'f99':
        if rv == None:
            rv=3.1
        aHaebv = f99_alebv(ha,rv)
        aHbebv = f99_alebv(hb,rv)
        return -2.5*np.log10(ha_hb/ha_hb_0)*(rv/(aHaebv-aHbebv))
# ---------------------------------------------------------------------------------------- 

def check():
    '''
    Make some quick plots to test the reddening functions of BRIAN, by comparing with
    the original papers. ANd make sure I did not mess things up ...
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
    '''
    x = np.arange(0.0001,7.5,0.1)
    lams = 1/x*1.e4
     
    plt.close(2)
    plt.figure(2, figsize=(10,8))
    gs = gridspec.GridSpec(2,1, height_ratios=[1,0.5], width_ratios=[1])
    gs.update(left=0.15,right=0.9,bottom=0.1,top=0.95,wspace=0.05,hspace=0.05 )
    
    r = 3.1
    f99_curve = f99_alebv(lams,r) - r
    ccm89_curve = (ccm89_alav(lams,r)*r) -r
    cal00_curve = cal00_ke(lams,4.05) - 4.05   
    fd05_curve = fd05_elvebv(lams)
    
    ax1 = plt.subplot(gs[0,0])
    ax1.plot(x,f99_curve,'r-', linewidth=1,label=r'F99, R$_V$='+np.str(r))
    ax1.plot(x,ccm89_curve,'g-',label=r'CCM89, R$_V$='+np.str(r))
    ax1.plot(x,cal00_curve,'b-',label=r'Cal00, R$_V$=4.05')
    ax1.plot(x,fd05_curve,'k--',linewidth=2, label='FD05, R$_V$$^{A}$=4.5')
    ax1.grid(True)
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'E($\lambda$-V)/E(B-V)')
    ax1.legend(loc='lower right',fontsize=15)
    ax1.set_xlim([0.5,4.5])
    ax1.set_ylim([-5,5])
    
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
       
    plt.show()
    '''
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
    #fd05_curve = fd05_elvebv(lams)+4.5
    
    ax1 = plt.subplot(gs[0,0])
    ax1.plot(np.log10(lams/1.e4),f99_2,'g:', linewidth=1,label=r'F99, R$_V$=2.0')
    ax1.plot(np.log10(lams/1.e4),f99_31,'g-', linewidth=1,label=r'F99, R$_V$=3.1')
    ax1.plot(np.log10(lams/1.e4),f99_5,'g--', linewidth=1,label=r'F99, R$_V$=5.0')
    ax1.plot(np.log10(lams/1.e4),ccm89_2,'r:', linewidth=1,label=r'CCM89, R$_V$=2.0')
    ax1.plot(np.log10(lams/1.e4),ccm89_31,'r-', linewidth=1,label=r'CCM89, R$_V$=3.1')
    ax1.plot(np.log10(lams/1.e4),ccm89_5,'r--', linewidth=1,label=r'CCM89, R$_V$=5.0')
    ax1.plot(np.log10(lams/1.e4),cal00_curve,'b-',label=r'Cal00, R$_V$=4.05')
    #ax1.plot(np.log10(lams/1.e4),fd05_curve,'k--',linewidth=2, label='FD05, R$_{V}$$^{A}$=4.5')

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


