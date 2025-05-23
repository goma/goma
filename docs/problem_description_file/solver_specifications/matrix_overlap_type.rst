***********************
Matrix Overlap Type
***********************

::

	Matrix overlap type = {standard | symmetric}

-----------------------
Description / Usage
-----------------------

This card selects the kind of matrix overlap that occurs (for parallel computations).
Valid options are:

standard
    The local processor considers only its own estimate for any unknown;
    results from adjacent processors are ignored. This is the default.
symmetric
    The local processor adds its own estimate together with estimates from
    adjacent processors, retaining symmetry of preconditioners if a symmetric
    technique is being employed.

If the *Matrix ovelap type* card is omitted, the default is **standard.**

------------
Examples
------------

Following is a sample card:
::

	Matrix overlap type= symmetric

-------------------------
Technical Discussion
-------------------------

This optional card determines how overlapping subdomain solver results are combined
when different processors derive different estimates for the same solution unknown.

This overlap option is moot for serial problems whose data decomposition is trivial.



