***************************
**Electrical Permittivity**
***************************

::

   Electrical Permittivity = {model_name} {float} []

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for electrical permittivity. There is
currently one option, so {model_name} must be either **CONSTANT**. Definitions of the
input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for constant electrical permittivity.                                                     |
|                 |                                                                                                            |
|                 | * <float> -the value of electrical permittivity.                                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+

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