***************************
Magnetic Permeability
***************************

::

   Magnetic Permeability = {model_name} {float} []

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for electrical permittivity.
Definitions of the input parameters are as follows:

CONSTANT
   Name of the model for constant electrical permittivity.                                                     |
      * <float> -the value of electrical permittivity.                                                           |
COMPLEX_CONSTANT
   complex permittivity (for Time-Harmonic Maxwell with Nedelec elements)
      * <float1> real part
      * <float2> imaginary part

------------
**Examples**
------------

Following are sample cards:

::

   Magnetic Permeability = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

This card is utilized to set the magnetic permeability for time-harmonic Maxwell problems



--------------
**References**
--------------

No References.
