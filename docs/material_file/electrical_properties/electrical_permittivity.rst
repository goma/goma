***************************
Electrical Permittivity
***************************

::

   Electrical Permittivity = {model_name} {float} []

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
REFRACTIVE_INDEX
   Calculate permittivity using the refractive index and extinction index cards (for Time-Harmonic Maxwell with Nedelec elements)
   .. math::

      \epsilon_r = (n - ik)^2

   Where :math:`n` is the refractive index and :math:`k` is the extinction index

------------
**Examples**
------------

Following are sample cards:

::

   Electrical Permittivity = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

This card is utilized to set the electrical permittivity for electrostatic problems.



--------------
**References**
--------------

No References.