# -*- coding: utf-8 -*-
#
# This file contains tools for the BRIAN routines to create pretty plots seamlessly.
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import aplpy

# ---------------------------------------------------------------------------------------- 
     
import matplotlib as mpl
# Use tex - a pain to setup, but it looks great when it works !
mpl.rc('MacOSX')
mpl.rc('font',**{'family':'sans-serif', 'serif':['Computer Modern Serif'], 
                 'sans-serif':['Helvetica'], 'size':20, 
                 'weight':500, 'variant':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':12, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm', 
                     'bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+\
                             r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'}) # dashed | solid
# ----------------------------------------------------------------------------------------      

def make_2Dplot(fn, # path to the data (complete!)
                ext = 1, # Which extension am I looking for ?
                ofn = '', # Savefig filename
                contours = False, # Draw any contours 
                stretch = 'arcsinh',
                vmin = None, 
                vmax = None,
                cmap = None,
                ):
    '''
    This function generates an image from a 2-D fits image. 
    '''
    
    # Plot the BW image using the 2d image to indicate the projection
    plt.close(1)
    fig1 = plt.figure(1, figsize=(10,9))
    
    ax1 = aplpy.FITSFigure(fn, hdu=ext, figure=fig1, north=True)
    
    if cmap:
        ax1.show_colorscale(stretch=stretch, vmax = vmax, vmin=vmin, 
                            cmap=cmap)
    
    else:
        ax1.show_grayscale(invert=True, stretch=stretch,vmax = vmax, vmin=vmin)

    ax1.set_tick_color('k')

    ax1.add_grid()
 
    if contours:
        ax1.show_contour(fn, hdu=ext, levels=contours,colors=['darkorange'],smooth=3) 
        #Smooth~seeing

    # Make it look pretty
    ax1.set_axis_labels(ylabel='Dec. (J2000)')
    ax1.set_axis_labels(xlabel='R.A. (J2000)')
    ax1.tick_labels.set_xformat('hh:mm:ss')
    ax1.tick_labels.set_yformat('dd:mm:ss')

    ax1.grid.set_color('k')
    ax1.grid.set_linestyle('dotted')

    # Save it
    fig1.savefig(ofn, bbox_inches='tight')
    plt.close(1)
    return True
# ----------------------------------------------------------------------------------------      
