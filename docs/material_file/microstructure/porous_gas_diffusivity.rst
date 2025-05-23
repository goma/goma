**************************
**Porous Gas Diffusivity**
**************************

::

   Porous Gas Diffusivity = {model_name} <integer> <float_list> [L2/t]

-----------------------
**Description / Usage**
-----------------------

This card sets the model for the porous gas diffusivity, or the diffusion coefficient for
diffusive species flux in the gas phase of a porous medium. It is applicable to media
types **POROUS_UNSATURATED** and **POROUS_TWO_PHASE** (see *Media Type*
card).

Definitions of the input options for {model_name} and the <integer> and <float>
parameters fro each model are as follows:

+----------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT** <integer>      |the name for the constant diffusivity model.                                         |
|<float1>                    |                                                                                     |
|                            | * <integer> - phase/component; always set to zero until a multicomponent capability |
|                            |   exists                                                                            |
|                            | * <float1> - D, Diffusivity [L2/t]                                                  |
+----------------------------+-------------------------------------------------------------------------------------+
|**POROUS** <integer>        |the name for a microstructure dependent porous medium model.                         |
|<float1> <float2> <float3>  |                                                                                     |
|<float4> <float5>           | * <integer> - phase/component; always set to zero until a multicomponent capability |
|                            |   exists                                                                            |
|                            | * <float1> - Dvo, binary diffusion coefficient in free space [L2/t]                 |
|                            | * <float2> - τ, tortuosity of the matrix skeleton                                   |
|                            | * <float3> - P*gas, reference gas phase pressure                                    |
|                            | * <float4> - T0, reference temperature.                                             |
|                            | * <float5> - n, exponent on the temperature dependence (see below).                 |
+----------------------------+-------------------------------------------------------------------------------------+

For two-phase or unsaturated flow in a porous medium, the diffusivity calculated by
this model is the diffusivity of solvent vapor through the gas phase in the pore-space
(see Martinez, 1995).

------------
**Examples**
------------

::

   Porous Gas Diffusivity = POROUS 0 1.e-5 0.5 1.e+6 25.0 3

See the equation below for the diffusivity model that this card represents.

-------------------------
**Technical Discussion**
-------------------------

The generalized flux of liquid phase solvent, in both gas and liquid phases, contains a
term that accounts for diffusion of the liquid solvent species as gas vapor (see
references below). That flux is as follows:

.. figure:: /figures/418_goma_physics.png
	:align: center
	:width: 90%

If the media type is **POROUS_TWO_PHASE**, this expression is divided by

.. figure:: /figures/419_goma_physics.png
	:align: center
	:width: 90%

and if in addition it is temperature dependent, this expression is multiplied by

.. figure:: /figures/420_goma_physics.png
	:align: center
	:width: 90%



--------------
**References**
--------------

GT-008.2: Porous Media Capabilities/Tutorial for GOMA. User Guidance for Saturated
Porous Penetration Problems, August 11, 1999, P. R. Schunk

GT-009.3: GOMA’s Capabilities for Partially Saturated Flow in Porous Media,
September 1, 2002, P. R. Schunk

SAND94-0379: “Formulation and Numerical Analysis of Nonisothermal Multiphase
Flow in Porous Media”, Sandia Technical Report, Martinez, M. J., 1995

SAND96-2149: Drying in Deformable Partially-Saturated Porous Media: Sol-Gel
Coatings, Cairncross, R. A., P. R. Schunk, K. S. Chen, S. S. Prakash, J. Samuel, A. J.
Hurd and C. Brinker (September 1996)