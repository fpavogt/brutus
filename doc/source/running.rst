
Running brian
===================

The spirit of brian is that each user can choose, depending on the object at hand & the
quality of the data, what processing steps are warranted. These are governed by the
``brian_execute.py`` file. In parallel, the user can define all relevant parameters using
the ``brian_params.py`` file. 

For each project, I suggest making a local copy of both files, while leaving an original 
copy at ``brian/exec_scripts/`` for safekeeping. 


Setting parameters with ``brian_params.py``
--------------------------------------------

All parameters (with a scientific role) relevant to the brian functions can be specified
inside ``brian_param.py``, via the dictionary ``brian_params``.

.. todo::
   
    A detailed list of what everything is/does would be nice.


Once everything is set as it should be, run the script ::

    >>> run brian_params.py

This exports all the parameters to the pickle file ``pickle_fn``, defined in 
``brian_params.py``.  


Launching the post-processing with ``brian_execute.py``
--------------------------------------------------------

brian is designed to batch process all the individual spectra inside a given IFU 
individually, exploiting multiple cpus (when available) to gain speed. 

The different steps to be executed in sequence by brian are defined in 
``brian_execute.py``, via the ``proc_steps`` list of dictionaries. 


Steps can be toogle on/off via the ``run`` keyword, which can be set to ``True`` or 
``False``. Some steps also accept specific parameters, such as whether to save plots or 
not. One should note that unlike the parameters defined in ``brian_params.py``, the 
parameters present in ``brian_execute.py`` will NOT influence the final scientific output 
itself, only the processing speed and/or the type of output, etc ...

Once all the steps are setup and enabled as required, you are ready to start crunching 
numbers::

    >>> run brian_execute pickle_fn
    
where ``pickle_fn`` is the filename defined in ``brian_params.py`` for exporting the 
various parameters. 

If all goes according to plan, brian will start processing your data, one step after the 
other. **Some steps will take time (up to several hours) !** Just be patient - progress 
indicators will help you monitor each task, if you set ``'verbose': True`` inside 
``brian_params.py``. 


