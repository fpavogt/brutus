.. _changelog:

Changelog
=========

v0.1.1, April 2016, F.P.A. Vogt
 - replaced '0.1' with '__version__' when saving fits headers
 - separated the continuum fitting from the construction of the continuum data cube, for
   clarity and consistency with the emission lines
 - added local copy of mpfit as brian_mpfit.py for portability
 - refined plot legends and colorbar ticks
 - started catching numpy warnings in known locations, triggered by all-nan's spaxels
 - implemented a dictionary of created filenames (in brian_metadata), to better keep track 
   of what is created, with much less hard-coded filenames
 - bugfixes
 
v0.1.0, April 2016, F.P.A. Vogt
 - created overall code structure inspired by pywifes
 - implemented continuum fitting using LOWESS approach
 - implemented modular emission line fitting via mpfit
 - added safety measures to catch KeyboardInterrupt during multiprocessing tasks.
 - created preliminary documentation structure 
 - initial upload to Github
 - added license info
 - added semi-automatic detection of structures, and ability to interact to add/remove
   apertures of varying radii
 - added creation of "aperture spectra cube" based on a given aperture list