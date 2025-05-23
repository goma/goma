**************************
Moment Weight Function
**************************

::

   Moment Weight Function = {GALERKIN | SUPG } <float>

-----------------------
Description / Usage
-----------------------

This card specifies the weight function to be used on the weighted residual of the
Moment equations. Definitions of the input
parameters are as follows:

GALERKIN     
    Name of the model for the weight functions for a full Galerkin formulation. This is the default when this 
    card is absent.
    * <float> - set but unused for Galerkin
SUPG
    Name of the model for the weight functions for a streamwise upwinded Petrov-Galerkin formulation.

    * <float> - the value of the weight function, a number between 0. and 1.; a value of 1. corresponds to a
      full **SUPG**

------------
Examples
------------

The following is a sample input card:

::

   Moment Weight Function = GALERKIN 0.0

   Moment Weight Function = SUPG 1.0

-------------------------
Technical Discussion
-------------------------


--------------
References
--------------

No References.