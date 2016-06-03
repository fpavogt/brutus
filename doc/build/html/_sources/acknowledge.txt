
Acknowledging brutus
====================

Only use lower case letters when mentioning brutus, and always include the release number, 
e.g.:

    brutus |release|  

brutus will be described in detail in 

    **Vogt et al.**, in prep.

If you find brutus useful for your research, please cite this reference accordingly. 

brutus uses several packages that **should also be acknowledged in their own right.** 
The following Tex-formatted acknowledgment is one way to do so::

    This research has made use of \textsc{brutus}, a Python module to process data cubes 
    from integral field spectrographs hosted at \url{http://fpavogt.github.io/brutus/}. 
    \textsc{brutus} relies on \textsc{statsmodel} (Seabold & Perktold 2010), 
    \textsc{ppxf} (Cappellari & Emsellem 2004), \textsc{fit_kinematic_pa} as described in 
    Appendix C of Krajnovic et al. (2006), \textsc{matplotlib} (Hunter 2007), 
    \textsc{astropy}, a community-developed core Python package for Astronomy 
    (Astropy Collaboration et al., 2013), \textsc{photutils}, an affiliated package of 
    \textsc{astropy} for photometry, \textsc{aplpy}, an open-source plotting package for 
    Python hosted at \url{http://aplpy.github.com}, \textsc{montage}, funded by the 
    National Science Foundation under Grant Number ACI-1440620 and previously funded by 
    the National Aeronautics and Space Administration’s Earth Science Technology Office, 
    Computation Technologies Project, under Cooperative Agreement Number NCC5-626 between 
    NASA and the California Institute of Technology, and \textsc{mpfit}, a Python script 
    that uses the Levenberg-Marquardt technique (Moré 1978) to solve least-squares 
    problems, based on an original Fortran code part of the \textsc{minpack}-1 package.

Finally, you also ought to cite the following works, depending on your use of brutus:

    1) Cleveland(1979); 
        the reference for the Locally Weighted Scatterplot Smoothing (LOWESS) algorithm used 
        by brutus (via statsmodels) to fit the continuum.
    2) The default SSP models shipped with brutus are::
    
        Stellar population synthesis model predictions using the MILES stellar libraries 
        (Sánchez-Blázquez et al. 2006, Falcón-Barroso et al. 2011) at FWHM=2.5 Å, based on 
        the code presented in Vazdekis et al. (2010), using a Kroupa revised IMF with 
        slope 1.3 (Kroupa 2001), Padova +00  (Girardi et al. 2000) isochrones, metallicities 
        ranging from -2.32 to +0.22, and ages of 0.0631 to 17.7828 Gyr.
            
    3) The reddening laws:
        Either the Cardelli, Clayton & Mathis (1989) law, the Calzetti et al. (2000) law or 
        the theoretical model of a turbulent dust screen of **Fischera & Dopita (2005) 
        [default]** for the extragalactic attenuation corrections, and
        the **Fitzpatrick (1999) law [default]** for the galactic extinction.
        
        If you use the extinction values :math:`A_V` and :math:`A_B` from the NASA 
        Extragalactic Database (NED) to correct for the Galactic extinction [default], then
        to be thorough, you should mention that::
        
            The Galactic extinction is derived using NED from the Schlafly & Finkbeiner 
            (2011) recalibration of the Schlegel, Finkbeiner & Davis (1998) infrared-based 
            dust map. The map is based on dust emission from COBE/DIRBE and IRAS/ISSA; 
            the recalibration assumes a Fitzpatrick (1999) reddening law with Rv = 3.1 and 
            different source spectrum than Schlegel, Finkbeiner & Davis (1998).
        
        and you also ought to acknowledge NED itself::
        
            This research has made use of the NASA/IPAC Extragalactic Database (NED) 
            which is operated by the Jet Propulsion Laboratory, California Institute of 
            Technology, under contract with the National Aeronautics and Space Administration. 
        
References:
 - `Astropy Collaboration et al. (2013) <http://cdsads.u-strasbg.fr/abs/2013A%26A...558A..33A>`_
 - `Cardelli, Clayton & Mathis (1989) <http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_
 - `Calzetti et al. (2000) <http://adsabs.harvard.edu/abs/2000ApJ...533..682C>`_
 - `Cappellari & Emsellem (2004) <http://cdsads.u-strasbg.fr/abs/2004PASP..116..138C>`_
 - Cleveland (1979)::
    
    @article{doi:10.1080/01621459.1979.10481038,
        author = { William S.   Cleveland },
        title = {Robust Locally Weighted Regression and Smoothing Scatterplots},
        journal = {Journal of the American Statistical Association},
        volume = {74},
        number = {368},
        pages = {829-836},
        year = {1979},
        doi = {10.1080/01621459.1979.10481038},
    }
 
 - `Falcón-Barroso et al. (2011) <http://adsabs.harvard.edu/abs/2011A%26A...532A..95F>`_
 - `Fischera & Dopita (2005) <http://adsabs.harvard.edu/abs/2005ApJ...619..340F>`_
 - `Girardi et al. (2000) <http://adsabs.harvard.edu/abs/2000A%26AS..141..371G>`_
 - `Hunter (2007) <http://cdsads.u-strasbg.fr/abs/2007CSE.....9...90H>`_    
 - `Krajnovic et al. (2006) <http://adsabs.harvard.edu/abs/2006MNRAS.366..787K>`_
 - `Kroupa 2001 <http://adsabs.harvard.edu/abs/2001MNRAS.322..231K>`_
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
    
 - `Sánchez-Blázquez et al. (2006) <http://adsabs.harvard.edu/abs/2006MNRAS.371..703S>`_
 - Seabold & Perktold (2010)::
 
    @inproceedings{seabold2010,
        title={Statsmodels: Econometric and statistical modeling with python},
        author={Seabold, Skipper and Perktold, Josef},
        booktitle={9th Python in Science Conference},
        year={2010},
    }
 - `Schlafly & Finkbeiner (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737..103S>`_  
 - `Schlegel, Finkbeiner & Davis (1998) <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_ 
 - `Vazdekis et al. (2010) <http://adsabs.harvard.edu/doi/10.1111/j.1365-2966.2010.16407.x>`_
