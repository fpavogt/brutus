
Acknowledging brian
====================

Only use lower case letters when mentioning brian, and always include the release number, 
e.g.:

    brian |release|  

brian will be described in detail in 

    **Vogt et al.**, in prep.

If you find brian useful for your research, please cite this reference accordingly. 

brian uses several packages that **should also be acknowledged in their own right.** 
The following Tex-formatted acknowledgment is one way to do so::

    This research has made use of \textsc{brian}, a Python module to process data cubes 
    from integral field spectrographs. \textsc{brian} relies on \textsc{statsmodel} 
    (Seabold & Perktold 2010), \textsc{ppxf} (Cappellari & Emsellem 2004), 
    \textsc{matplotlib} (Hunter 2007), \textsc{astropy}, a community-developed core Python 
    package for Astronomy (Astropy Collaboration et al., 2013), \textsc{photutils}, an 
    affiliated package of \textsc{astropy} for photometry, \textsc{aplpy}, an open-source 
    plotting package for Python hosted at \url{http://aplpy.github.com}, \textsc{montage}, 
    funded by the National Science Foundation under Grant Number ACI-1440620 and 
    previously funded by the National Aeronautics and Space Administration’s Earth Science 
    Technology Office, Computation Technologies Project, under Cooperative Agreement 
    Number NCC5-626 between NASA and the California Institute of Technology, and 
    \textsc{mpfit}, a Python script that uses the Levenberg-Marquardt technique 
    (Moré 1978) to solve least-squares problems, based on an original Fortran code part of 
    the \textsc{minpack}-1 package.

References:
 - `Astropy Collaboration et al. (2013) <http://cdsads.u-strasbg.fr/abs/2013A%26A...558A..33A>`_
 - `Cappellari & Emsellem (2004) <http://cdsads.u-strasbg.fr/abs/2004PASP..116..138C>`_
 - `Hunter (2007) <http://cdsads.u-strasbg.fr/abs/2007CSE.....9...90H>`_    
 - Seabold & Perktold (2010)::
 
    @inproceedings{seabold2010,
        title={Statsmodels: Econometric and statistical modeling with python},
        author={Seabold, Skipper and Perktold, Josef},
        booktitle={9th Python in Science Conference},
        year={2010},
    }
 - Moré (1978)::
 
    @inbook{more1978,
        author={Mor{\'e}, Jorge J.},
        editor={Watson, G. A.},
        chapter={The Levenberg-Marquardt algorithm: Implementation and theory},
        title={Numerical Analysis: Proceedings of the Biennial Conference 
               Held at Dundee, June 28--July 1, 1977},
        year={1978},
        publisher={Springer Berlin Heidelberg},
        address={Berlin, Heidelberg    },
        pages={105--116},
        isbn={978-3-540-35972-2},
        doi={10.1007/BFb0067700},
        url={http://dx.doi.org/10.1007/BFb0067700}
    }

