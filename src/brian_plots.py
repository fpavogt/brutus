# -*- coding: utf-8 -*-
#
# This file contains tools for the BRIAN routines to create pretty plots seamlessly.
#
# Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
# ----------------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec 
import aplpy
import photutils
from astropy.stats import sigma_clipped_stats

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
                cblabel = None,
                cbticks=None,
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
    
    ax1.add_colorbar()
    ax1.colorbar.set_width(0.2)
    ax1.colorbar.set_location('top')
    if cblabel:
        ax1.colorbar.set_axis_label_text(cblabel)
        ax1.colorbar.set_axis_label_pad(10)
    if cbticks:
        ax1.colorbar.set_ticks(cbticks)
    ax1.grid.set_color('k')
    ax1.grid.set_linestyle('dotted')

    # Save it
    fig1.savefig(ofn, bbox_inches='tight')
    plt.close(1)
    return True
# ----------------------------------------------------------------------------------------      

class ApManager(object):
    '''
    This class is designed to add and remove apertures of different radii to a
    matplotlib plot interactively. It handles right-clicks, left-clicks and keyPressEvents
    I understand half-of-it, the rest being some black magic from the web.
    '''
    def __init__(self, fig, ax, points, rs):
        self.axes = ax
        self.fig = fig
        self.canvas = ax.figure.canvas
        self.points = points
        self.xs = list(points.get_xdata())
        self.ys = list(points.get_ydata())
        self.rs = list(rs)
         
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.keyPress = self.canvas.mpl_connect('key_press_event', self.onKeyPress)
        self.pick=self.canvas.mpl_connect('pick_event', self.onpick)
         
        self.radKey = 3. # Default aperture radius for new apertures, in pixels.
        self.pickEvent = False
    # -------------------------------------------------------------------------------------      

    def onpress(self, event):
        '''
        Define what happens when a click happens on the matplotlib window.
        '''
        #print 'Press Event received: %s' % event
        if self.pickEvent: # Are we also having a pickEvent ? Then, let it deal with it.
            self.pickEvent = False
            return
        
        # Black magic - no idea what this is. removed for now.
        #if self.canvas.widgetlock.locked(): return 
        
        # Did the click happen outside the axes ?
        if event.inaxes is None: return
            
        # Is anything selected in the toolbar (e.g. zoom mode ?) 
        # Only proceed if not.
        tb = plt.get_current_fig_manager().toolbar
        if event.inaxes!=self.points.axes: 
            return
            
        # Left click ? Then proceed to add a new aperture.
        if tb.mode == '' and event.button == 1:
            # First, add the aperture center, which I use for the picking.
            # Round the result to an integer.
            self.xs.append(np.round(event.xdata))
            self.ys.append(np.round(event.ydata))
            self.rs.append(self.radKey)
            self.points.set_data(self.xs, self.ys)
            #self.points.figure.canvas.draw()
            
            # Then, the aperture footprint, which I do NOT pick, and simply update it
            # on the way. It's a patch added to the ax.   
            my_ap_scatter(self.axes, [np.round(event.xdata)], [np.round(event.ydata)], 
                       [self.radKey], 
                       facecolor='none',edgecolor='tomato')
             
            # And update the plot           
            self.fig.canvas.draw()
    # ------------------------------------------------------------------------------------     

    def onKeyPress(self, event):
        '''
        Allows to change the aperture size by typing a number. Currently in pixels.
        Number needs to be typed while the cursor is on the plot window.
        '''
        # Just to be nice, in case the user missed a pushed button on the toolbar
        # Test if everything is free, issue a warning otherwise.
        tb = plt.get_current_fig_manager().toolbar
        if tb.mode != '':
            sys.stdout.write('\r   Warning: the plot toolbar is busy            ')
            sys.stdout.flush()   
            return
            
        # Using the numeric keys causes trouble with the zoom level. Instead, use "u" and
        # "d" for incremental changes to the aperture size !
        if event.key == 'u': 
            self.radKey += 1.
            sys.stdout.write('\r   New aperture radius: %i pixels ' % self.radKey)
            sys.stdout.flush()
            
        if event.key == 'd': 
            if self.radKey > 1:
                self.radKey -= 1.
            sys.stdout.write('\r   New aperture radius: %i pixels ' % self.radKey)
            sys.stdout.flush()
    # ------------------------------------------------------------------------------------     

    def onpick(self, event):
        '''
        Defines what happens with pickEvent. Namely, remove the near-by aperture.
        '''
        #print 'Pick Event received: %s' % event
        self.pickEvent = True # Used to avoid clashes between onpick and onpress. I think.
        # If right click, then remove the aperture. 
        if event.mouseevent.button == 3:
            index = event.ind
            # Remove the xs, ys and rs value.
            self.xs.pop(index)
            self.ys.pop(index)
            self.rs.pop(index)
            
            # Remove the peak point
            self.points.set_data(self.xs,self.ys)
            #self.points.figure.canvas.draw()
            
            # Also remove the aperture
            self.axes.patches.pop(index)
            
            # Update the figure
            self.fig.canvas.draw()
# ----------------------------------------------------------------------------------------      
            
def ap_outline(x0,y0,radius):
    '''
    Construct the exact aperture outline, for all pixels within "radius" of a given pixel.
    This works well for R <5, but starts to look weird afterwards. Basically, my function
    starts missing pixels, and the area start looking like a diamond. The alternative is
    to define all the otuline polygons by hand - but I really do not feel like doing this.
    Is there a way to construct them with a smarter routine than below ?
    This function is left here for legacy purposes and remind me of what I did. But it
    is not used by brian at the moment.
    '''
    nu = np.linspace(-radius-0.5, +radius+0.5,2*(radius+1))
    nd = np.linspace(+radius+0.5, -radius-0.5,2*(radius+1))
    
    xs = np.append([x for pair in zip(nu,nu) for x in pair][1:-1],
                [x for pair in zip(nd,nd) for x in pair][1:-1])
    
    ys = np.append( xs[-len(xs)/4:],xs[:-len(xs)/4])
    
    return zip(np.array(xs)+x0,np.array(ys)+y0)
# ----------------------------------------------------------------------------------------      
    
def my_ap_scatter(ax, xs, ys, rs, **kwargs):
    '''
    This function is a "pseudo-plt.scatter" function, but that instead use the plt.Circle
    routine to show individual apertures of different radii. Importantly, the radii
    are thus expressed in data units, and not in points !
    '''
                            
    for (x, y, r) in zip(xs, ys, rs):
        circle = plt.Circle((x,y), radius=r, **kwargs)
        # TODO: Instead of a Circle, could I show the outline of which pixels will be 
        # included or not ? The function below works, but not for radius >5, when it 
        # looks weird (points start missing).
        #circle = plt.Polygon(ap_outline(x,y,r), **kwargs)
        ax.add_patch(circle)
        
    return True
# ----------------------------------------------------------------------------------------      
        
def build_ap_list(data, 
                  start_aps = None,
                  radius = 3.0
                  ):
    '''
    This function finds local maxima in an image, and then offers the user the ability
    to check/modify the found peaks manually. It associate a circular aperture to each 
    peak, used later to to extract the integrated spectra of these areas.  
    '''
    
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5) 
    
    # I tried to play with daofind at some point, with some limited success. 
    # It is very sensitive to varying background, and was not ideal to find HII regions
    # at the core of galaxies. 
    # The next lines is left here for possible future usage, and remind me that I tried.
    # sources = daofind(data, fwhm=5.0, threshold=1.*std) 

    # For now, these parameters are fixed. Until I see whether playing with them at the 
    # user level can help or not.
    threshold =  (1.0 * std) 
    #myf = np.array([[0,1,1,1,0],
    #                [1,1,1,1,1],
    #                [1,1,1,1,1],
    #                [1,1,1,1,1],
    #                [0,1,1,1,0]])               
    
    if not(start_aps):
        # Find the local peak in the data.
        tbl = photutils.find_peaks(data, threshold, box_size=5, subpixel=False)
        xs = tbl['x_peak']
        ys = tbl['y_peak']
        rs = np.zeros(len(xs))+radius
    else:
        # User provided us with apertures already. No need to start from scratch
        xs,ys,rs = zip(*start_aps)
    
    # Start the plotting
    plt.close(90)
    fig = plt.figure(90, figsize=(10,8))
    gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
  
    ax = plt.subplot(gs[0,0])
 
    # Note, I will here be working in pixel space for simplicity. No need to bother
    # with WCS coefficient for now.
    ax.imshow(data, cmap='Greys', origin='lower', 
              norm=LogNorm(), vmin=1,
              interpolation='nearest')
              
    # Plot the peaks, and make them pickable !        
    line, = ax.plot(xs, ys, ls='none', color='tomato',
             marker='+', ms=10, lw=1.5, picker=4)
             
    # Plot the associated apertures, but DON'T make them pickable.
    # Construct the size of the apertures associated with each peak, using the 
    # default value from the user.
    # TODO: be smarter, and find the best aperture size for each peak. 
    # E.G. 2D gaussian fitting ?
    # Also for now, just use a Circle            
    my_ap_scatter(ax, xs, ys, rs, facecolor='none',edgecolor='tomato')  
    
    # Now, the 1-line black magic command that deals with the user interaction          
    apmanager = ApManager(fig, ax, line, rs)          

    ax.grid(True)
    ax.set_xlim((-10,np.shape(data)[1]+10))
    ax.set_ylim((-10,np.shape(data)[0]+10))
    plt.show()

    print '   Starting interactive mode:'
    print '   Aperture radius: 3 pixels'
    print '      [left-click]      : add an aperture' 
    print '      [right-click]     : remove an aperture'
    print '      [u]+cursor on plot: larger aperture '
    print '      [d]+cursor on plot: smaller aperture '                    
    
    # A little trick to wait for the user to be done before proceeding
    while True:
    
        letter = raw_input('   [m] to save and continue, [zzz] to NOT save and crash. \n')
        if letter in ['m','zzz']:
            break
        else:
            print '   [%s] unrecognized !' % letter
        
    plt.close()
    
    if letter =='m':
        return zip(apmanager.xs, apmanager.ys, apmanager.rs)
    else:
        return False
# ----------------------------------------------------------------------------------------      
    
    
    
    
    