*****************
**Charge Number**
*****************

::

   Charge Number = CONSTANT <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card is required when charged species are involved, e.g. when using the
**FICKIAN_CHARGED** or the **STEFAN_MAXWELL_CHARGED** *Diffusion
Constitutive Equation* card. It specifies the charge number (e.g., the charge number for
Ni2+ is 2, and that for SO2-is -2) of a species.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Model for specifying constant charge on species.                                     |
|                          |                                                                                     |
|                          | * <integer> - species number                                                        |
|                          | * <float> - charge number of the species                                            |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Sample usage for this card is shown belo

::

   Charge Number = CONSTANT 0 1.0

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.