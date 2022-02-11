********************
**Volume Expansion**
********************

::

   Volume Expansion = CONSTANT <float> [1/T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the coefficient of volume expansion in the
energy equation. This property is required for the **BOUSS** option on the *Navier-Stokes
Source* card. Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for a constant volume-expansion coefficient.                                              |
|                 |                                                                                                            |
|                 | * <float> - the value of the volume expansion coefficient.                                                 |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Volume Expansion = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

Warning: Please be careful that the Species Volume Expansion card is set
appropriately. If the **BOUSS** or **BOUSSINESQ** model is turned on on the Navier-
Stokes Source card, then both thermal and species volume expansion effects are
included if the coefficients are nonzero. .



--------------
**References**
--------------

No References.