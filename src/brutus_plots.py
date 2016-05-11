# -*- coding: utf-8 -*-
'''
 brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2016,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains tools for the brutus routines to create pretty plots seamlessly.

 Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec 
import aplpy
import photutils
from astropy.stats import sigma_clipped_stats

from brutus_metadata import alligator

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

def make_galred_plot(lams, alams, etau, Ab, Av, ofn):
    ''' Plots the extinction Alambda and flux correction factor for galactic extinction.
    
    :Args:
        lams: 1-D numpy array
              The wavelength nodes.
        alams: 1-D numpy array
               The corresponding values of Alambda.
        etau: 1-D numpy array
              The flux correction factor, i.e. F_(unextinct)/F_(observed).
        Ab: float
            The value of Ab.
        Av: float
            The value of Av.
    
    :Returns: True
              Always.
    '''
    
    plt.close(1)
    plt.figure(1,figsize=(10,7))
    gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
  
    ax = plt.subplot(gs[0,0])
    
    ax.plot(lams,alams,'k-',drawstyle='steps-mid')
    
    ax.grid(True)
    ax.set_xlabel(r'Wavelength [\AA]')
    ax.set_ylabel(r'A$_\lambda$ [mag]')
    
    ax2 = ax.twinx()
    ax2.plot(lams,etau,'-',c='firebrick',drawstyle='steps-mid')
    ax2.set_ylabel(r'F$_{\lambda,0}$/F$_{\lambda,obs}= e^{\tau}$', color='firebrick',
                   labelpad=10)
    
    ax2.tick_params(axis='y',colors='firebrick')
    ax.set_title(r'Galactic extinction correction')
    
    
    plt.savefig(ofn, bbox_inches='tight')
    plt.close() 
    


# ----------------------------------------------------------------------------------------      

def make_2Dplot(fn, # path to the data (complete!)
                ext = 1, # Which extension am I looking for ?
                ofn = 'plot.pdf', # Savefig filename
                contours = False, # Draw any contours 
                stretch = 'arcsinh',
                vmin = None, 
                vmax = None,
                cmap = None,
                cblabel = None,
                cbticks = None,
                ):
    '''Creates an image from a 2-D fits image with WCS info.
    
    :Args:
        fn: string
            The filename (+path!) fo the fits file to display.
        ext: int (default:1)
            The extension of the image to display. The data in this extension MUST be 2-D.
        ofn: string [default:'plot.pdf']
            The output filename (+path!)
        contours: bool,list of floats [default: False]
            If set, shows contours at said levels, e.g. [1,3,5,10.5]
        stretch: string [default:'arcsinh']
            The stretch of the image, fed to aplpy.FITSFigure. 'linear', 'arcsinh', ...
        vmin: float [default: None]
            If set, the lower bound of the colorbar
        vmax: float [default: None]
            If set, the upper bound of the colorbar
        cmap: string
            If set, the colormap to use. Use 'alligator' for the special brutus cmap.
        cblabel: string [default:None]
            If set, the label of the colorbar.
        cbticks: list of int [default: None]
            If set, the colorbar ticks.
            
    :Returns:
        out: True
             Always.
        
    :Notes:
        This function absolutely requires WCS coordinates, and a 2-D array.
    '''
    
    # Plot the BW image using the 2d image to indicate the projection
    plt.close(1)
    fig1 = plt.figure(1, figsize=(10,9))
    
    ax1 = aplpy.FITSFigure(fn, hdu=ext, figure=fig1, north=True)
    
    if cmap == 'alligator':
        ax1.show_colorscale(stretch=stretch, vmax = vmax, vmin=vmin, 
                            cmap=alligator)
    elif cmap:
        ax1.show_colorscale(stretch=stretch, vmax = vmax, vmin=vmin, 
                            cmap=cmap)
    else:
        ax1.show_grayscale(invert=True, stretch=stretch,vmax = vmax, vmin=vmin)

    ax1.set_tick_color('k')

    ax1.add_grid()
 
    if contours:
        ax1.show_contour(fn, hdu=ext, levels=contours,colors=['darkorange']) 
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
    ax1.set_nan_color((0.5,0.5,0.5))

    # Save it
    fig1.savefig(ofn, bbox_inches='tight')
    plt.close(1)
    return True
# ----------------------------------------------------------------------------------------      

class ApManager(object):
    ''' The class managing the interactive aperture selection plot.
    
    This class is designed to add and remove apertures of different radii to a
    matplotlib plot interactively. It handles right-clicks, left-clicks and keyPressEvents
    
    :Args: 
        fig: Figure instance from matplotlib
             The Figure we are working with.
        ax: Axes instance from matplotlib
            The axes we are working with.
        points: list instance from plt.plot
            The points to pick.
        rs: list of int
            The radius of the apertures, in pixels.
    
    :Notes:
        I understand 90% of this class, the rest being some black magic from the web.
    '''
    def __init__(self, fig, ax, points, rs):
        ''' Defines the specific elements required by the class.'''
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
        '''Defines what happens when an 'event' happens on the matplotlib window.
        
        This function allows to add an aperture of a given radius, when doing a left-click
        on a matplotlib plot window.
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
        '''Defines what happens when a key is pressed, while on a matplotlib window.
        
        This function allows to change the aperture size by typing 'u' (radius += 1 pixel)
        or 'd' (radius -= 1 pixel). Minimum size set to 1 pixel.
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
        '''Defines what happens with pickEvent occurs in a matplotlib window.
        
        This function removes a given aperture if the user right-clicks on it.
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
    ''' A function to construct the aperture outline as a function of the radius.
    
    :Args:
        x0: int
            x-coordinate of the aperture center, in pixels
        y0: int
            y-coordinate of the aperture center, in pixels
        radius: int
            aperture radius, in pixels.
    
    :Returns:
        out: 2-D numpy array
            The list of [xs,ys] coordinates of the aperture outline nodes, in pixels.
    
    :Notes:
       This function is in principle designed to construct the exact aperture outline for 
       all pixels within "radius" of a given pixel. This works well for R <5, but breaks 
       afterwards. Basically, the function starts missing pixels, and the area start 
       looking like a diamond. The alternative is to define all the otuline polygons by 
       hand - but I really do not feel like doing this manually. Is there a way to 
       construct them with a smarter routine than below ? This function is left here for 
       legacy purposes and remind me of what I did. But it is not used by brutus at the 
       moment. Apertures are all shown as perfect circles for now.
    '''
    nu = np.linspace(-radius-0.5, +radius+0.5,2*(radius+1))
    nd = np.linspace(+radius+0.5, -radius-0.5,2*(radius+1))
    
    xs = np.append([x for pair in zip(nu,nu) for x in pair][1:-1],
                [x for pair in zip(nd,nd) for x in pair][1:-1])
    
    ys = np.append( xs[-len(xs)/4:],xs[:-len(xs)/4])
    
    return zip(np.array(xs)+x0,np.array(ys)+y0)
# ----------------------------------------------------------------------------------------      
    
def my_ap_scatter(ax, xs, ys, rs, **kwargs):
    ''' A pseudo-plt.scatter function for easily displaying many apertures.
    
    This function use the plt.Circle routine to show individual apertures of different 
    radii. Importantly, the radii are thus expressed in data units, and not in points.
    
    :Args:
        ax: Axes instance from matplotlib
            The Figure axes to which to add the aperture.
        xs: list of int
            The x-coordinates of the aperture centers, in pixels.
        ys: list of int
            The y-coordinates of the aperture centers, in pixels.
        rs: list of int
            The radii of the apertures, in pixels.
    
    :Returns:
        out: True
             Always.
        
    :Notes:
        In principle, each aperture could be shown "exactly" using plt.Polygon and the 
        brutus_plots.ap_outline function. But the latter is not suitable for large 
        apertures. Until fixed, each apertuyre is shown as a perfect circle using 
        plt.Circle. 
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
                  radius = 3.0,
                  automatic_mode=True,
                  interactive_mode=False,
                  lam = None,
                  save_plot = None,
                  ):
    ''' Detects local maxima in an image, and allows the manual inspection of the result.
    
    This function finds local maxima in 2-D array, and then offers the user the ability
    to check/modify the found peaks manually. It associate a fixed circular aperture to 
    each peak (individually modifiable manually), used later to extract the integrated s
    pectra of these areas. 
    
    :Args:
        data: 2-D numpy array
            The raw data from which to find maxima.
        start_aps: list of list [default: None]
            A pre-existing list of apertures.
        radius: int [default: 3.0]
            The default radius of the apertures associated with the local maxima detected
            automatically.
        automatic_mode: bool [default: True]
                        Whether to detect peaks and assign apertures automatically or not.
        interactive_mode: bool [default:True]
                          Whether to inspect the apertures manually or not.
        lam: float [default: False]
             The wavelength of the data, to be added to the plot title if set.
        save_plot: string [default: None]
                   If set, the name of the file used to save the final aperture selection.
        
    :Returns:
        out: 2-D numpy array
             The array of apertures with shape (n,3), where n is the total number of
             apertures (1 per row), defined by (x,y,r) its center and radius (in pixels). 
    :Notes:
        If start_aps are provided, then automatic_mode is False.     
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
    
    if not(start_aps) and automatic_mode:
        # Find the local peak in the data.
        tbl = photutils.find_peaks(data, threshold, box_size=5, subpixel=False)
        xs = tbl['x_peak']
        ys = tbl['y_peak']
        rs = np.zeros(len(xs))+radius
    elif start_aps:
        # User provided us with apertures already. No need to start from scratch
        xs,ys,rs = zip(*start_aps)
    else:
        xs = np.zeros(0)
        ys = np.zeros(0)
        rs = np.ones(0)
        
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
    if interactive_mode:          
        apmanager = ApManager(fig, ax, line, rs)          

    ax.grid(True)
    ax.set_xlim((-10,np.shape(data)[1]+10))
    ax.set_xlabel(r'x [pixels]', labelpad = 10)
    ax.set_ylim((-10,np.shape(data)[0]+10))
    ax.set_ylabel(r'y [pixels]', labelpad = 10)
    ax.set_title(r'$\lambda$:%.2f' % lam)
    plt.show()

    if interactive_mode:
        print '   Starting interactive mode:'
        print '   Aperture radius: 3 pixels'
        print '      [left-click]      : add an aperture' 
        print '      [right-click]     : remove an aperture'
        print '      [u]+cursor on plot: larger aperture '
        print '      [d]+cursor on plot: smaller aperture '                    
    
        # A little trick to wait for the user to be done before proceeding
        while True:
    
            letter = raw_input('   [m] to save and continue, [zzz] to NOT save and '+
                               'crash (maybe). \n')
            if letter in ['m','zzz']:
                break
            else:
                print '   [%s] unrecognized !' % letter
        
        line.remove()
        plt.savefig(save_plot+np.str(lam)+'.pdf',bbox_inches='tight')            
        plt.close()
        
        if letter =='m':
            return zip(apmanager.xs, apmanager.ys, apmanager.rs)
        else:
            return False
    else:
        line.remove()
        plt.savefig(save_plot+np.str(lam)+'.pdf',bbox_inches='tight')            
        plt.close()
        return zip(xs,ys,rs)
# ----------------------------------------------------------------------------------------      

class SpecManager(object):
    '''The class managing the interactive inspection of the fitted data.
    
    This class is designed to select spaxels, and plot the associated fit.
    
    :Notes:
        This class is an evolution of brutus_plots.ApManager(). It will make anyone with
        experience in events handling cry of despair. It works for now, but it could be 
        faster, and with no doubt A LOT more elegant. Suggestions welcome.
    '''
    def __init__(self, fig, ax1a, ax1b, ax3a, ax3b, ax3c, ax3d,
                       patches, spectra, lams, data, lowess, ppxf, elines):

        self.fig = fig
        self.ax1a = ax1a
        self.ax1b = ax1b
        self.ax3a = ax3a
        self.ax3b = ax3b
        self.ax3c = ax3c
        self.ax3d = ax3d

        self.canvas = ax1a.figure.canvas
        
        self.patch1a = patches[0]
        self.patch1b = patches[1]
        
        self.spec3 = spectra[0]
        self.lowess3 = spectra[1]
        self.ppxf3 = spectra[2]
        self.fit3 = spectra[3]
        
        self.dspec3b = spectra[4]
        self.dspec3c = spectra[5]
        self.dspec3d = spectra[6]
        
        self.lams = lams
        self.data = data
        self.lowess = lowess
        self.ppxf = ppxf
        self.elines = elines
         
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def onpress(self, event):
        '''
        #Define what happens when a click happens on the matplotlib window.
        '''
        # Did the click happen outside the axes ?
        if event.inaxes is None: return
            
        # Is anything selected in the toolbar (e.g. zoom mode ?) 
        # Only proceed if not.
        tb = plt.get_current_fig_manager().toolbar

        # Left click on the images ? Then load a new spectrum from the chosen spaxel
        if tb.mode == '' and event.button == 1 and event.inaxes in [self.ax1a, self.ax1b]:
            # First, update the spaxel marker location
                
            self.patch1a.remove()
            symb = plt.Rectangle([np.int(event.xdata)-0.5,np.int(event.ydata)-0.5],
                                     1,1,angle=0,color='tomato')
            self.patch1a = self.ax1a.add_patch(symb)
            # ---
            self.patch1b.remove()
            symb = plt.Rectangle([np.int(event.xdata)-0.5,np.int(event.ydata)-0.5],
                                     1,1,angle=0,color='tomato')
            self.patch1b = self.ax1b.add_patch(symb)
            
            # Update the label
            self.ax1a.set_title('Spaxel: %s - %s' % (str(np.int(event.xdata)).zfill(3),
                                                    str(np.int(event.ydata)).zfill(3)),
                                                    fontsize=20)
                
            # Now, update the spectrum in the master spec plot
            for dat in [self.spec3, self.lowess3, self.ppxf3, self.fit3, self.dspec3b,
                        self.dspec3c, self.dspec3d]:
                dat[0].remove() # Erase   
            
            self.spec3 = self.ax3a.plot(self.lams,
                                       self.data[:,int(event.ydata),int(event.xdata)],
                                       'k-',drawstyle='steps-mid') # Add new data
                      
            
            self.ppxf3 = self.ax3a.plot(self.lams,
                                       self.ppxf[:,int(event.ydata),int(event.xdata)],
                                       '-',c='firebrick',drawstyle='steps-mid', lw=2)  
            self.lowess3 = self.ax3a.plot(self.lams,
                                       self.lowess[:,int(event.ydata),int(event.xdata)],
                                       '-',c='darkorange',drawstyle='steps-mid', lw=2) 
            self.fit3 = self.ax3a.plot(self.lams,
                                       self.lowess[:,int(event.ydata),int(event.xdata)]+
                                       self.elines[:,int(event.ydata),int(event.xdata)],
                                       '-',c='royalblue',drawstyle='steps-mid', lw=1) 
            
            self.dspec3b = self.ax3b.plot(self.lams,
                                       self.data[:,int(event.ydata),int(event.xdata)]-
                                       self.lowess[:,int(event.ydata),int(event.xdata)],
                                       '-',c='k',drawstyle='steps-mid', lw=1)     
            self.dspec3c = self.ax3c.plot(self.lams,
                                       self.data[:,int(event.ydata),int(event.xdata)]-
                                       self.ppxf[:,int(event.ydata),int(event.xdata)],
                                       '-',c='k',drawstyle='steps-mid', lw=1) 
            self.dspec3d = self.ax3d.plot(self.lams,
                                       self.data[:,int(event.ydata),int(event.xdata)]-
                                       self.lowess[:,int(event.ydata),int(event.xdata)]-
                                       self.elines[:,int(event.ydata),int(event.xdata)],
                                       '-',c='k',drawstyle='steps-mid', lw=1)                                             
            
            # And reset the axes y-range automatically
            self.ax3a.relim()
            self.ax3a.autoscale_view()
             
            # And update the plot 
            for ax in [self.ax3a,self.ax3b, self.ax3c, self.ax3d, self.ax1a, self.ax1b]:          
                ax.figure.canvas.draw()
            
        # Left click on the dspec area ? Double the axis range
        elif tb.mode == '' and event.button == 1 and event.inaxes in [self.ax3b, 
                                                                      self.ax3c,
                                                                      self.ax3d]:
                                                                    
            currlim = self.ax3b.get_ylim()   
            for ax in [self.ax3b,self.ax3c,self.ax3d]: 
                ax.set_ylim((currlim[0]*2,currlim[1]*2))
                ax.set_yticks((currlim[0],currlim[1]))
                ax.figure.canvas.draw()   
                                                            
        # Right click on the dspec area ? Half the axis range
        elif tb.mode == '' and event.button == 3 and event.inaxes in [self.ax3b, 
                                                                      self.ax3c,
                                                                      self.ax3d]:
                                                                    
            currlim = self.ax3b.get_ylim()   
            for ax in [self.ax3b,self.ax3c,self.ax3d]: 
                ax.set_ylim((currlim[0]/2.,currlim[1]/2.))
                ax.set_yticks((currlim[0]/4.,currlim[1]/4.))
                ax.figure.canvas.draw()   
                                           
            
# ----------------------------------------------------------------------------------------      

def inspect_spaxels(lams, data,lowess, ppxf,elines,map,vmap, irange, vrange, ofn=False):
    ''' Interactive inspection fo the brutus fit output.
    
    This function generates an interactive plot allowing to select individual spaxels 
    and visualize how well/bad they were fitted.
    
    :Args:
        lams: 1-D array
              The wavelength array of the data.
        data: 3-D numpy array
            The raw data.
        lowess: 3-D numpy array
            The fitted lowess continuum cube.
        ppxf: 3-D numpy array
            The fitted ppxf continuum cube.
        elines: 3-D numpy array
            The pure emission line cube.
        map: 2-D numpy array
            The line inetnsity map.
        vmap: 2-D numpy array
            The gas velocity map.
        irange: list of int
                The range of the colorbar for the intensity plot.
        vrange: list of int       
                The range of the colorbar for the velocity dispersion plot.
        ofn: string [default: None]       
              The filename (+path!) of the plot saved.
              
    :Returns:
        out: True
             Always.  
    '''
    # Starting location for the plot. Middle of the cube
    
    nx = int(np.shape(data[0])[1]/2.)
    ny = int(np.shape(data[0])[0]/2.)
     
    plt.close(91)
    fig = plt.figure(91, figsize=(16,9))

    # Here, create multiple instance of gridspec, to have "total" control.
    # First for two maps, to show the location of the spaxel.
    gs1 = gridspec.GridSpec(2,1, height_ratios=[1,1], width_ratios=[1])
    gs1.update(left=0.63, right=0.98, bottom=0.1, top=0.93, wspace=0, hspace=0.05)

    # Then 4 spectra plots, to show the various fits.
    gs3 = gridspec.GridSpec(4,1, height_ratios=[1,0.5,0.5,0.5], width_ratios=[1])
    gs3.update(left=0.1, right=0.63, bottom=0.1, top=0.96, wspace=0, hspace=0.05)
    
    # Start with the 2-D maps
    ax1b = plt.subplot(gs1[1,0])    
    ax1a = plt.subplot(gs1[0,0],sharex=ax1b, sharey=ax1b) # Here, sharex/y = zoom is tied!
    ax1a.set_adjustable('box-forced')
    ax1b.set_adjustable('box-forced')
     
    # Here, I will need to allow for an automatic range set for the colorbars 
    im1 = ax1a.imshow(map, cmap=alligator, origin='lower', norm=LogNorm(), vmin=irange[0],
                      vmax=irange[1], 
                      interpolation='nearest')         
    im2 = ax1b.imshow(vmap, cmap=alligator, origin='lower', vmin=vrange[0],vmax=vrange[1], 
                      interpolation='nearest')               
     
    # A red dot for showin which spaxel we're at. 
    symb = plt.Rectangle([nx-0.5,ny-0.5],1,1,angle=0,color='r')
    patch1a = ax1a.add_patch(symb)
    symb = plt.Rectangle([nx-0.5,ny-0.5],1,1,angle=0,color='r')
    patch1b = ax1b.add_patch(symb)
    patches = [patch1a,patch1b]
    
    # Also write it down, for clarity
    ax1a.set_title('Spaxel: %s - %s' % (str(np.int(nx)).zfill(3), 
                                        str(np.int(ny)).zfill(3)),
                                        fontsize=20) 
    # Deal with the axis labels, etc...           
    ax1a.set_xticklabels([])
    for ax in [ax1a, ax1b]:
        ax.yaxis.tick_right()
        ax.yaxis.label_position='right'  # Put this on the right.
        ax.set_ylabel(r'y [pixels]', labelpad=30, fontsize=20) 
        if ax in [ax1b]:
            ax.set_xlabel(r'x [pixels]', fontsize=20)
            
    # Now, deal with the various spectra
    ax3d = plt.subplot(gs3[3,0])
    ax3a = plt.subplot(gs3[0,0], sharex=ax3d) # Tie the zoom, here again.
    ax3b = plt.subplot(gs3[1,0], sharex=ax3d) # Tie the zoom, here again.
    ax3c = plt.subplot(gs3[2,0], sharex=ax3d) # Tie the zoom, here again.
    
    # Hide the axis ticks hen not necessary
    for ax in [ax3a,ax3b,ax3c]:
        plt.setp(ax.get_xticklabels(), visible=False)
    
    # Start the actual plotting
    spec3 = ax3a.plot(lams,data[:,ny,nx],'k-',drawstyle='steps-mid',lw=2)
    ppxf3 = ax3a.plot(lams,ppxf[:,ny,nx],'-',c='firebrick', drawstyle='steps-mid', lw=2)   
    lowess3 = ax3a.plot(lams,lowess[:,ny,nx],'-',c='darkorange', drawstyle='steps-mid',
                                                                                   lw=2)  
    fit3 = ax3a.plot(lams,lowess[:,ny,nx]+elines[:,ny,nx],'-',c='royalblue', 
                                                                   drawstyle='steps-mid') 
                                                                                                                           
    dspec3b = ax3b.plot(lams,data[:,ny,nx]-lowess[:,ny,nx],'k-',drawstyle='steps-mid')    
    dspec3c = ax3c.plot(lams,data[:,ny,nx]-ppxf[:,ny,nx],'k-',drawstyle='steps-mid')    
    dspec3d = ax3d.plot(lams,data[:,ny,nx]-lowess[:,ny,nx]-elines[:,ny,nx],'k-',
                                                                   drawstyle='steps-mid')    
    spectra = [spec3,ppxf3,lowess3,fit3,dspec3b,dspec3c,dspec3d]

    # Make things look pretty
    for ax in [ax1a,ax1b,ax3a,ax3b,ax3c,ax3d]:
        ax.grid(True)
        if ax in [ax3a,ax3b,ax3c,ax3d]:
            ax.set_xlim((np.min(lams),np.max(lams)))
        if ax in [ax3b,ax3c,ax3d]:
            ax.set_ylim((-20,20))
            ax.set_yticks([-10,10])
            ax.axhline(0,c='r')
        if ax in [ax3d]:
            ax.set_xlabel(r'Observed wavelength [\AA]')
    # Add labels
    ax3a.set_ylabel(r'F$_\lambda$ [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]', fontsize=15)
    ax3b.set_ylabel(r'F$_\lambda$ - lowess ', fontsize=15)
    ax3c.set_ylabel(r'F$_\lambda$ - ppxf ', fontsize=15)
    ax3d.set_ylabel(r'F$_\lambda$ - fit ', fontsize=15)
    
    # Launch the interactive black magic        
    specmanager = SpecManager(fig,ax1a,ax1b,ax3a,ax3b,ax3c,ax3d,
                              patches, spectra, lams,data, lowess, ppxf, elines)       
    plt.show()
    
    # A little trick to wait for the user to be done before proceeding
    while True:
    
        letter = raw_input('  [v]+[some-text] to save the current view, [z] to end.\n')
        if letter[0] == 'v':
            plt.savefig(ofn+letter+'.pdf',bbox_inches='tight')
        elif letter == 'z':
            break
        else:
            print '   [%s] unrecognized !' % letter
                 
    plt.close(91)
        
    return True