**************************
Moment Shock Function
**************************

::

   Moment Shock Function = {NONE | YZBETA } <float>

-----------------------
Description / Usage
-----------------------

This card specifies the shock capturing function to be used on the 
Moment equations. Definitions of the input
parameters are as follows:

NONE     
    No shock capturing is applied. This is the default when this 
    card is absent
YZBETA
    A YZBETA shock capturing term is applied.

    * <float> - the value of the weight function, a number between 0. and 1.; a value of 1. corresponds to a
      full YZBETA shock capturing

------------
Examples
------------

The following is a sample input card:

::

   Moment Weight Function = NONE (default)

   Moment Weight Function = YZBETA 1.0

-------------------------
Technical Discussion
-------------------------


--------------
References
--------------

No References.