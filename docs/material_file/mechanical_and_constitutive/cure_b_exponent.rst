***************
Cure B Exponent
***************

::

   Cure B Exponent = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the B exponent of the **CURE**, **EPOXY**, and
**FILLED_EPOXY** model options of the *Liquid Constitutive Equation* card.
Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the *B* exponent.                                                                     |
|                 |                                                                                                            |
|                 | * <float> - the value of *B*.                                                                              |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the cure *B* exponent to 0.1:

::

   Cure B Exponent = CONSTANT 0.1.

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



