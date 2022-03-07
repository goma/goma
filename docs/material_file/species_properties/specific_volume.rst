*******************
**Specific Volume**
*******************

::

   Specific Volume = CONSTANT <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card specifies the specific volume of a species. It is required when polymersolvent
vapor-liquid phase equilibrium is considered at the material boundaries.
Specific volume is used to convert units of mass fraction to volume fraction in species
material balance.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Specific volume model type.                                                          |
|                          |                                                                                     |
|                          | * <integer> - species number                                                        |
|                          | * <float> - pure component specific volume                                          |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Specific Volume = CONSTANT 0 1.154734411

-------------------------
**Technical Discussion**
-------------------------

This is the place where pure component density (inverse of specific volume)
information is entered in the material property. When Flory-Huggins vapor-liquid
equilibrium model was first developed in *Goma*, the equations were based on volume
fractions, not mass fractions. In order to convert these units, the specific volume
parameter is required for each component in the mixture.

This card is used only in conjunction with Flory-Huggins nonideal liquid activity
model for polymer-solvent mixtures. This occurs when two types of BCs are specified:
1) when *VL_POLY* is specified at an discontinuous internal boundary and 2) when
FLORY model under *YFLUX_EQUIL* boundary card is specified.



--------------
**References**
--------------

No References.