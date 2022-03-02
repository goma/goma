*************************
**Polymer Time Constant**
*************************

::

   Polymer Time Constant = {model_name} <float> [t]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the polymer time constant associated with the *Polymer
Constitutive Equation* card. It is a required card for the **OLDROYDB, GIESEKUS**
and **PTT** options. Definitions of the input parameters are as follows:

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

If the polymer time constant varies with properties, it must do so in the same way as the
polymer viscosity; thus, the model on this card must be the same as the model selected
on the *Polymer Viscosity* card.

All three models are described in detail in Bird, et. al. (1987).

------------
**Examples**
------------

The following is a sample card that sets the polymer time constant to 1.0:

::

   Polymer Time Constant = CONSTANT 1.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

Bird, R. B., Armstrong, R. C., and Hassager, O. 1987. Dynamics of Polymeric Liquids,
2nd ed., Wiley, New York, Vol. 1.