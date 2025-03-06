***********************
Species Shock Capturing
***********************

::

   Species Shock Capturing = {model_name} <float>

-------------------
Description / Usage
-------------------

This optional card is used to specify the weight functions to be used on the weighted
residual of the species convective diffusion equations. For high Peclet number cases,
you may want to use a Petrov-Galerkin formulation rather than a Galerkin formulation.

{model_name}       
    Name of the formulation model. Valid entries are:

    * NONE
    * MIXED

<float> 
    independent scaling to apply, 1 is full shock capturing and 0 is no shock capturing.

When this card is absent, the default {model_name} is **NONE**.

--------
Examples
--------

The following is a sample input card:

::

   Species Shock Capturing = MIXED 0.5

--------------------
Technical Discussion
--------------------

No Discussion.



----------
References
----------

No References.
