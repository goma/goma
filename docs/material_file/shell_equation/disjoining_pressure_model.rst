*****************************
**Disjoining Pressure Model**
*****************************

::

   Disjoining Pressure Model = {model_name} <floatlist>

-----------------------
**Description / Usage**
-----------------------

This card takes the specification of the disjoining pressure model for the film-flow
equation capability, specifically the shell_filmh equation. This function
specifies the disjoining pressure in the unit of force per area. Currently four models for
{model_name} are permissible:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |This model specifies a constant disjoining pressure. This option requires one        |
|                          |floating point values                                                                |
|                          |                                                                                     |
|                          | * <float1> is the evaporation rate in the unit of length/time                       |
+--------------------------+-------------------------------------------------------------------------------------+

**ONE_TERM** This model specifies disjoining pressure and the input parameters
(floats). This model only employs the repulsion part of the van der Waals force. The
functional form is:

.. figure:: /figures/495_goma_physics.png
	:align: center
	:width: 90%

* <float1> is the equilibrium liquid-solid contact angle θe

* <float2> is exponent n and it should satisfy n >1

* <float3> is the precursor film thickness h*

**TWO_TERM** This model specifies disjoining pressure and the input parameters
(floats). Here, the model only employs both repulsion and attraction part of the van der
Waals force. The functional form is:

.. figure:: /figures/496_goma_physics.png
	:align: center
	:width: 90%

where

.. figure:: /figures/497_goma_physics.png
	:align: center
	:width: 90%

* <float1> is the equilibrium liquid-solid contact angle 0e

* <float2> is exponent n corresponding to the repulsive part of the
  van der Waals force. It should satisfy n >1

* <float3> is exponent m corresponding to the attractive part of the
  van der Waals force. It should satisfy m > n since
  the attractive part acts in longer range than the repulsive
  one.

* <float4> is the precursor film thickness h

* <float5> is the parameter a describing relative importance of the
  attractive part to the repulsive part. Typically, a is
  chosen to be 0 <α <1 in order to achieve
  more numerical stability.

**TWO_TERM_EXT_CA** This model is identical with **TWO_TERM** except that it
uses contact angle from an external field identifies as THETA.

------------
**Examples**
------------

Following is a sample card:

::

   Disjoining Pressure Model = TWO_TERM 120.3 2 1.0e-4 0.1

This results in disjoining pressure with contact angle of 120, repulsive exponent of 3,
attractive repulsion of 2, precursor film thickness of 1.0e-4, and relative importance of
attractive part of 0.1.

-------------------------
**Technical Discussion**
-------------------------

A thorough discussion of disjoining pressure can be found in Teletzke et al (1987). The
premultiplying constant is related to contact angle and surface tension by
balancing capillary and disjoining force where the wetting line meets the precursor film.
See Schwartz (1998) for further detail.



--------------
**References**
--------------

Leonard W. Schwartz, R. Valery Roy, Richard R. Eley, and Stanislaw Petrash,
“Dewetting Patterns in a Drying Liquid Film”, Journal of Colloid and Interface
Science 234, 363–374 (2001)

Teletzke, G. F., Davis, H. T., and Scriven, L. E., “How liquids spread on solids”, Chem.
Eng. Comm., 55, pp 41-81 (1987).