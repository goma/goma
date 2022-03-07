********************************
**Non-volatile Specific Volume**
********************************

::

   Non-volatile Specific Volume = CONSTANT <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card specifies the specific volume of a species when the species is implicit in the
mixture. This means that in most problems involving n+1 species, only n species are
independent; i.e.,

.. figure:: /figures/461_goma_physics.png
	:align: center
	:width: 90%

It is required when Flory-Huggins vapor-liquid phase equilibrium is considered at the
material boundaries, as used in *VL_POLY* and in **FLORY** under *FLUX_EQUIL*. This
is used to convert units of mass fraction to mole fraction in species material balance.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Model for converting mass to mole fractions.                                         |
|                          |                                                                                     |
|                          | * <integer> - species number                                                        |
|                          | * <float> - specific volume of the non-volatile species, usually the n+1 component  |
|                          |   in *Goma* convention.                                                             |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is an example card:

::

   Non-volatile Specific Volume = CONSTANT 2 0.855e-3

This example shows that two species are solved in the *Goma* problem explicitly:
species 0 and species 1.

-------------------------
**Technical Discussion**
-------------------------

In the current set up, species balance in *Goma* considers the species to be independent
of each other. However, the mass or volume fractions of all species must add up to
unity in any mixtures. This means that some properties of the last species must be
entered in the material file although that component is not solved explicitly in the
problem. This is the case for molecular weight, molar volume, and specific volume
specifications, which are required for calculating Flory-Huggins liquid activity.



--------------
**References**
--------------

No References.