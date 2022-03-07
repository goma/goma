***********************
**Solidus Temperature**
***********************

::

   Solidus Temperature = CONSTANT <float> [T]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the solidus temperature. Definitions of the
input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for the solidus temperature.                                                              |
|                 |                                                                                                            |
|                 | * <float> - the value of the solidus, Ts .                                                                 |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card:

::

   Solidus Temperature = CONSTANT 1.

-------------------------
**Technical Discussion**
-------------------------

This card is required when using the **ENTHALPY** option on the *Heat Capacity* card.



--------------
**References**
--------------

No References.