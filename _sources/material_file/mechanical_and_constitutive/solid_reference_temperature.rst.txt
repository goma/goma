***************************
Solid Reference Temperature
***************************

::

   Solid Reference Temperature = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the solid reference temperature used in the
thermal expansion of solid materials. Definitions of the input parameters are as
follows:

+-----------------+-------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the reference temperature.                               |
+-----------------+-------------------------------------------------------------------------------+
|<float>          |A floating point number that is the value of the solid reference temperature,  |
|                 |:math:`T_{ref}` .                                                              |
+-----------------+-------------------------------------------------------------------------------+


------------
**Examples**
------------

The following is a sample card:

::

   Solid Reference Temperature = CONSTANT 90.0

-------------------------
**Technical Discussion**
-------------------------

See the *Solid Thermal Expansion* card for a discussion of the use of this property in the
linear thermoelasticity of solids.



--------------
**References**
--------------

No References.