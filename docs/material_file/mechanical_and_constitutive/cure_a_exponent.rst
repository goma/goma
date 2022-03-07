***************
Cure A Exponent
***************

::

   Cure A Exponent = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the A exponent of the **CURE**, **EPOXY**, and
**FILLED_EPOXY** model options of the *Liquid Constitutive Equation* card.
Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the *A* exponent.                                                                     |
|                 |                                                                                                            |
|                 | * <float> - the value of *A*.                                                                              |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the cure *A* exponent to 1.0:

::

   Cure A Exponent = CONSTANT 1.0

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



