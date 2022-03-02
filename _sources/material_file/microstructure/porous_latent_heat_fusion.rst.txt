*****************************
**Porous Latent Heat Fusion**
*****************************

::

   Porous Latent Heat Fusion = CONSTANT <integer> <float> [E/M]

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the model for the latent heat of fusion (or freezing)
for each liquid solvent species in a partially saturated porous media flow problem, viz.
*Media Type* card set to **POROUS_UNSATURATED** or **POROUS_TWO_PHASE**.
As of 6/13/2002, we only allow single liquid phase solvent and the porous enthalpy
equation is being tested. Definitions of the input parameters are as follows:

+----------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**                |Name of the constant latent heat of fusion model.                                    |
|                            |                                                                                     |
|                            | * <integer> - the species equation of liquid phase solvent; MUST BE SET TO ZERO for |
|                            |   now.                                                                              |
|                            | * <float> - the value of the latent heat of fusion.                                 |
+----------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Porous Latent Heat Fusion = CONSTANT 0 1000.2

See the equation below for the diffusivity model that this card represents.

-------------------------
**Technical Discussion**
-------------------------

First order phase change involves the adsorption or expulsion of heat. This thermal
effect is modeled through the porous energy equation (see *EQ* cards; this equation was
under development and testing as this manual was being assembled) with a source term
that depends on the evaporation/condensation rate. *Fusion* implies a liquid to solid
transition. It is envisioned that this card will someday be used for porous flow in mushy
zones of solidifying metals, or the freezing of water in a porous solid.



--------------
**References**
--------------

No References.