**************
Cure Gel Point
**************

::

   Cure Gel Point = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the :math:`α_g` parameter for the **CURE**, **EPOXY**,
and **FILLED_EPOXY** model options of the *Liquid Constitutive Equation* card.
Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the :math:`α_g` parameter.                                                            |
|                 |                                                                                                            |
|                 | * <float> - the value of :math:`α_g`, which is the extent of reaction at the gel point of a polymerizing   |
|                 |   system.                                                                                                  |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the cure gel point to 0.75:

::

   Cure Gel Point = CONSTANT 0.75

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



