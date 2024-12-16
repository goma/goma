**************************
Moment Time Integration
**************************

::

   Moment Time Integration = {STANDARD | TAYLOR_GALERKIN }

-----------------------
Description / Usage
-----------------------

This card specifies the time integration method to be used for time-stepping of the
Moment equations. Definitions of the input
parameters are as follows:

STANDARD     
    Name of the method for time-stepping in the standard formulation. This is the default when this card is absent.

TAYLOR_GALERKIN
    Name of the method for time-stepping in using the Implicit Taylor-Galerkin formulation.

------------
Examples
------------

The following is a sample input card:

::

   Moment Time Integration = STANDARD

   Moment Time Integration = TAYLOR_GALERKIN


-------------------------
Technical Discussion
-------------------------


--------------
References
--------------

No References.
