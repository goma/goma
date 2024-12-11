**************************
Moment Time Integration Method Function
**************************

::

   Moment Time Integration = {STANDARD | TAYLOR_GALERKIN } <float>

-----------------------
Description / Usage
-----------------------

This card specifies the time integration method to be used for time-stepping of the
Moment equations. Definitions of the input
parameters are as follows:

STANDARD     
    Name of the method for time-stepping in the standard formulation. This is the default when this 
    card is absent.
    * <float> - set but unused for Standard
TAYLOR_GALERKIN
    Name of the method for time-stepping in using the Implicit Taylor-Galerkin formulation.
    * <float> - a number between 0. and 1.; a value of 1. corresponds to **TAYLOR_GALERKIN**

------------
Examples
------------

The following is a sample input card:

::

   Moment Time Integration = STANDARD 0.0

   Moment Time Integration = TAYLOR_GALERKIN 1.0


-------------------------
Technical Discussion
-------------------------


--------------
References
--------------

No References.
