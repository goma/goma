*************************
**Reference Temperature**
*************************

::

   Reference Temperature = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the reference temperature, which is required
by the **BOUSS** option on the *Navier-Stokes Source* card and by the **BINGHAM** option
on the *Liquid Constitutive Equation* card. Definitions of the input parameters are as
follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for a constant reference temperature.                                                     |
|                 |                                                                                                            |
|                 | * <float> - the value of the reference temperature.                                                        |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Reference Temperature = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.