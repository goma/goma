********************
**Molecular Weight**
********************

::

   Molecular Weight = CONSTANT <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card specifies the molecular weight of a species. It is required when the Stefan-
Maxwell flux model is used in modeling multicomponent transport of neutral or
charged species. It is also required when vapor-liquid phase equilibrium is considered
at the material boundaries. Molecular weight is used to convert units of mass fraction to
mole fraction in a species material balance.

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |Molecular weight model type.                                                         |
|                          |                                                                                     |
|                          | * <integer> - species number                                                        |
|                          | * <float> - molecular weight of the species                                         |
+--------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Molecular Weight = CONSTANT 0 6.939

-------------------------
**Technical Discussion**
-------------------------

This card originated from the development of a multicomponent diffusion model based
on the Stefan-Maxwell equation. However, it has been generalized to include problems
where mole fractions are necessary for the consideration of phase equilibria. For
example, when *YFLUX_EQUIL* is invoked in the input deck, an equilibrium problem is
solved rigorously which requires gas and liquid mole fractions. The conversion from
mass fraction to mole fraction requires molecular weight information.



--------------
**References**
--------------

No References.