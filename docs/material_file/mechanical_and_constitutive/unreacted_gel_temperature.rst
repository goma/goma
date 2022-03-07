*************************
Unreacted Gel Temperature
*************************

::

   Unreacted Gel Temperature = CONSTANT <float>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the unreacted gel temperature parameter for
the **FILLED_EPOXY** fluid option of the *Liquid Constitutive Equation* card.

Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the unreacted gel temperature.                                                        |
|                 |                                                                                                            |
|                 | * <float> - the value of the unreacted gel temperature,:math:`T_{g0}`.                                     |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the unreacted gel temperature to 273.0:

::

   Power Law Exponent = CONSTANT 273.0

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.




