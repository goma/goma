*********************
**Polymer Viscosity**
*********************

::

   Polymer Viscosity = {model_name} <float> [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the polymer viscosity associated with the model set in the
*Polymer Constitutive Equation* card. This is a required card for the **OLDROYDB**,
**GIESEKUS** and **PTT** models.

Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|{model_name}     |Permissible names for the viscosity model are **CONSTANT, POWER_LAW** and **CARREAU**.                      |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |a simple constant viscosity, Newtonian fluid.                                                               |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**POWER_LAW**    |a power-law model                                                                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**CARREAU**      |a Carreau strain-rate thinning or thickening relation                                                       |
+-----------------+------------------------------------------------------------------------------------------------------------+

Input parameters are not identified for the latter two models as they have not worked
since the multimode port. They could be made to work again if the proper tweaking is
done to the input parser, but are not currently functional.

------------
**Examples**
------------

The following is a sample card that sets the polymer viscosity to 8000.0:

::

   Polymer Viscosity = CONSTANT 8000.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------
