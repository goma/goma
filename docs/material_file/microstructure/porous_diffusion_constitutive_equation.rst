******************************************
**Porous Diffusion Constitutive Equation**
******************************************

::

   Porous Diffusion Constitutive Equation = {model_name}

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the species diffusion model for the gas phase in a
porous medium. Just now there is only one option, but plans are to expand the options
to include multicomponent diffusion models (cf. *Diffusion Constitutive Equation* card).
It is important to note that this model specification only applies to the gas phase of each
component. Liquid phase species diffusive transport has not been implemented as of
12/19/01.

Definitions of the input parameters are as follows, with only a single permissible value:

+-------------------+-------------------------------------------------------------------------------------+
|**DARCY_FICKIAN**  |Name of the model for the diffusion constitutive equation in the porous gas phase.   |
+-------------------+-------------------------------------------------------------------------------------+

This model simply implies that gas species can be transported relative to the solid
skeleton phase not only by a pressure gradient, as in Darcy’s law, but also by Fickian
diffusion.

------------
**Examples**
------------

The following sample input card uses the APREPRO variable model_name (which is
set to **DARCY_FICKIAN**.

::

   Porous Diffusion Constitutive Equation = {model_name}

-------------------------
**Technical Discussion**
-------------------------

Currently, the **DARCY_FICKIAN** model is the only option for the porous diffusion
equation and it only applies to one phase. When this card is parsed, it is contained in a
solvent species loop. When we allow more than one volatile species, we will eventually
allow for other diffusion constitutive equation models, e.g. of the Stefan-Maxwell type.
Also, we will have to build a phase dependence into this card, as the diffusion law may
be different in the liquid and in the gas. Right now, we do not allow for diffusion
transport (viz. by a chemical potential or concentration gradient) in the liquid phase of
a porous medium. Please consult references below for theoretical discussion.



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)