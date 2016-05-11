
Installing brutus
===================

The most recent release of brutus is available for download from its `Github repository
<https://github.com/fpavogt/brutus/releases>`_. 

.. note::
    You can also get the current (unreleased and possibly buggy) code from 
    `here <https://github.com/fpavogt/brutus/archive/master.zip>`_.

Installing brutus merely requires to let your Python installation know about its existence. 
Specifically:

1. Unzip the compressed files (downloaded from the link above), 
2. place them anywhere you like, and 
3. add ``brutus/src/`` to your Python path. 

In my case (on MAC OSX), my ``.bash_profile`` would look like : 
::

      export PYTHONPATH=$PYTHONPATH:/Users/fvogt/Tools/Python/fpav_pylib/brutus/src/


Github users can fork the brutus repository if they want to get access to the 
latest updates not yet formally released. *Push requests* for bug fixes and new features 
are welcome and will be examined in detail. 

.. _requirements:
      
Requirements
------------
The basic packages below are required for brutus to work properly:

* numpy (1.9.0 or or above)
* scipy (0.14.0 or above)
* matplotlib (1.4.2 or above)
* `aplpy <https://aplpy.github.io/>`_
* montage (`download <http://montage.ipac.caltech.edu/docs/download.html>`_ | `install instructions <http://montage.ipac.caltech.edu/docs/build.html>`_)
* `statsmodel <http://statsmodels.sourceforge.net/>`_
* `astropy <http://www.astropy.org/>`_
* `photutils <http://photutils.readthedocs.io/en/latest/>`_

Testing the installation
------------------------

.. todo::
   Create an elegant little thing to test that brutus is working fine. For now, just try to
   ``import brutus``.

.. _troubleshooting:

Troubleshooting
---------------

If something goes bad, here is what you should do (mind the order !):

1. 
    Are all the required packages up-to-date ? Check the :ref:`requirements`.
  
2. 
    Still no luck ? Check the :ref:`faq`.
  
3. 
    Still no luck ? Check if the issue has been `reported before 
    <https://github.com/fpavogt/brutus/issues>`_.
  
4. 
    Still no luck ? Please `submit a new issue <https://github.com/fpavogt/brutus/issues>`_ 
    on the Github repository of the project.
 
    Provide as much detail as possible (error message, minimal example able to reproduce 
    the error, operating system, Python version, multiprocessing state, etc ...).

    .. note::
       Submitting a Github issue is the best way for you to get help rapidly, for everyone 
       to keep track of the problems that need solving, and for future users to see what 
       changes have been made over time (and look for existing solutions to their problem 
       which may be the same as yours). Submitting a new issue on Github is rapid and easy, 
       but if you are really against doing it (why would you ?), you can always email 
       frederic.vogt@alumni.anu.edu.au for help. 

 