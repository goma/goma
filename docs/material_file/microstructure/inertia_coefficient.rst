***********************
**Inertia Coefficient**
***********************

::

   Inertia Coefficient = CONSTANT <float> []

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the inertia coefficient cˆ in the Brinkman
formulation for flow through porous media, viz. see **POROUS_BRINKMAN** option
on *Media Type* card. Detailed discussion of this coefficient can be found by consulting
the references below. Definitions of the input parameters are as follows:

+-----------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**                 |Name of the model for the inertia coefficient.                                       |
|                             |                                                                                     |
|                             | * <float> - The value of the inertia coefficient.                                   |
+-----------------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample input card that produces a weighting coefficient of 1.0 on the
inertial term in the **POROUS_BRINKMAN** equations.

::

   Inertia Coefficient= CONSTANT 1.0

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
