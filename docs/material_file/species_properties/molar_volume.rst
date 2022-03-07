****************
**Molar Volume**
****************

::

   Molar Volume = CONSTANT <integer> <float> [L3/mole]

-----------------------
**Description / Usage**
-----------------------

This card is referred when molar based equilibrium models are used on the boundaries,
such as *VL_POLY*. The float value specified is necessary for converting mass fractions
to mole fractions.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Model for converting mass to mole fractions.                                         |
|                          |                                                                                     |
|                          | * <integer> - species number                                                        |
|                          | * <float> - molar volume of the species. [L3/mole]                                  |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

An example usage for this card:

::

   Molar Volume = CONSTANT 0 1.

-------------------------
**Technical Discussion**
-------------------------

The same conversion from mass fraction to mole fraction can be obtained through
specification of the *Molecular Weight* and *Specific Volume*. The redundancy, which will
be allowed to remain, arose through simultaneous additions to the code by developers
working on different projects.



--------------
**References**
--------------

No References.