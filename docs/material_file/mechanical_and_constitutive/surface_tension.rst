*******************
**Surface Tension**
*******************

::

   Surface Tension = {CONSTANT | DILATION | USER} <float_list> [M/t2]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the interfacial surface tension of the fluid which enters into
the *CAPILLARY* boundary condition and *CAP_ENDFORCE* boundary condition cards.
The surface tension, albeit a property of an interface and not of a bulk material, is
sometimes influenced by thermophysical phenomena associated with a material, hence
the inclusion of this card in the material file. It should be mentioned that this card is
optional, and if it does not appear the surface tension is taken off the aforementioned
boundary condition cards. *PLEASE see the important technical discussion below if you
plan on using this card*. Definitions of the input parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**CONSTANT**     |Name of the model for a constant value of interfacial surface tension of the fluid:                         |
|                 |                                                                                                            |
|                 | * <float1> - constant value of surface tension                                                             |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**DILATION**     |Name of a surface tension model that depends on mesh dilation (only useful if the free surface is           |
|                 |constrained to be a material surface both normally and tangentially, see Schwartz ,et. al. 1996). The model |
|                 |mathematically is                                                                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**USER**         |A user-defined surface tension model that is defined in the user-supplied routine usr_surface_tension in the|
|                 |file user_mp.c. This model will have an arbitrary number of user-defined parameters (<float1> to floatn>).  |
+-----------------+------------------------------------------------------------------------------------------------------------+

*WARNING: When specifying surface tension on this card, be sure the surface
tension (multiplier) on the boundary condition CAPILLARY card is set to 1. In other
words, the value of surface tension on the boundary condition cards is multiplied
with the value on this card before the calculation is carried out.*

------------
**Examples**
------------

Following is a sample card:

::

   Surface Tension = DILATION 70.0 1.

-------------------------
**Technical Discussion**
-------------------------

Please read and understand the warning issued above regarding the proper place to
specify surface tension. Basically, for constant surface tension models, it is a good idea
to leave this card out and simply enter the proper surface tension value for the current
surface on the boundary condition cards **CAPILLARY** and **CAP_ENDFORCE**. For
variable models, please set the surface tension values on these BC cards to 1.0, and
then handle your model through this card. The surface tension is a thermodynamic
property of the interface and actually depends on the chemical composition of the
fluids (or fluid/solids) of the bounding phases. The property controls the importance of
the capillary stress jump on a curved interface on the hydrodynamics of the flow and
the meniscus position and motion.



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche