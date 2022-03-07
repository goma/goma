****************************
**Flowing Liquid Viscosity**
****************************

::

   FlowingLiquid Viscosity = CONSTANT <float> [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the viscosity of liquid flowing through pores
with the Brinkman model of flow through porous media, viz. see *Media Type* card with
**POROUS_BRINKMAN** option. In the Brinkman model, the viscosity input through
the *Viscosity* card is used as the Brinkman viscosity, and the viscosity input through this
card is used in determining the hydraulic resistance. Detailed discussion of these two
viscosities can be found by consulting the references below.

Definitions of the input parameters are as follows:

+-----------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**                 |Name for the constant viscosity model.                                               |
|                             |                                                                                     |
|                             | * <float> - The value of the viscosity.                                             |
+-----------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:

::

   Flowing Liquid Viscosity = CONSTANT 101.0

This card is only applicable to the **POROUS_BRINKMAN** media type and results in a
hydraulic resistance viscosity of 101.0.

-------------------------
**Technical Discussion**
-------------------------

See references below for discussion on use of this card.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

Gartling, D. K., C. E. Hickox and R. C. Givler 1996. "Simulations of Coupled Viscous
and Porous Flow Problems", Comp. Fluid Dynamics, 7, 23-48.