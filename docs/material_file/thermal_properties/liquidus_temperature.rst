************************
**Liquidus Temperature**
************************

::

   Liquidus Temperature = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the liquidus temperature. Definitions of the
input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the liquidus temperature.                                                             |
|                 |                                                                                                            |
|                 | * <float> - the value of the liquidus, Tl .                                                                |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Liquidus Temperature = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

This card is required when using the **ENTHALPY** option on the *Heat Capacity* card.



--------------
**References**
--------------

No References.