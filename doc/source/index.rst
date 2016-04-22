.. brian documentation master file, created by
   sphinx-quickstart on Fri Apr 22 11:23:38 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

brian |release|
==================

.. warning::
   These pages describe brian |release|. The module is currently under construction, and 
   NOT yet ready for deployment. For more information, send an email to 
   frederic.vogt@alumni.anu.edu.au.

   You can track the latest changes on the Github repository of the code: 
   https://github.com/fpavogt/brian

   See the :ref:`changelog` for more details.

brian is a Python module designed to post-process datacubes from Integral Field Spectrographs,
and in particular MUSE on the VLT/UT-4 at Paranal Observatory. brian is designed as a 
comprehensive set of routines to fit the stellar continuum and emission lines efficiently,
and automatically, with no input from the user - other than an inital sets of parameters.

brian is build in a modular fashion, with a structure inspired by 
`pywifes <http://pywifes.github.io/pipeline/>`_, the Python data reduction module for the 
WiFeS integral field spectrograph 
(`Childress+, 2014 <http://adsabs.harvard.edu/abs/2014Ap%26SS.349..617C>`_). Namely, brian 
contains a series of core modules, called by via a modular "execution sequence" defined by 
the user, alongside dedicated parameters. 

In particular, brian is designed to allow users to easily insert extra processing steps,
or by-pass others at will. brian uses the multiprocessing module, and can exploit up to 
300 cpus at once. brian also exploits existing tools with proven track-records developed 
by the community and individuals, including `astropy <http://www.astropy.org/>`_, 
`pyneb <http://www.iac.es/proyecto/PyNeb/>`_ and 
`ppxf <http://www-astro.physics.ox.ac.uk/~mxc/software/>`_. 

Contents:
------------------
.. toctree::
   :maxdepth: 1

   Home <self>
   installation
   running
   faq
   changelog
   acknowledge
   modules/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


----

Copyright notice:
 
This file is part of the brian Python module.
The brian Python module is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
 
The brian Python module is distributed in the hope that it will be useful, but without any 
warranty; without even the implied warranty of merchantability or fitness for a particular 
purpose.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with the brian 
Python module.  If not, see http://www.gnu.org/licenses/ .
 
