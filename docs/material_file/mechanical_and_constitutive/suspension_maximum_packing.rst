**************************
Suspension Maximum Packing
**************************

::

   Suspension Maximum Packing = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the :math:`C{max}` parameter of the **SUSPENSION**
and **FILLED_EPOXY** model options of the *Liquid Constitutive Equation* card.
Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for suspension maximum packing.                                                           |
|                 |                                                                                                            |
|                 | * <float> - the value of :math:`C{max}`, which is the mass fraction at which the suspension begins to act  |
|                 |   as a solid.                                                                                              |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the suspension maximum packing:

::

   Suspension Maximum Packing = CONSTANT 0.68.

-------------------------
**Technical Discussion**
-------------------------

See Description/Usage for *Liquid Constitutive Equation* card.



