****************
**PF_CAPILLARY**
****************

::

	BC = PF_CAPILLARY LS <integer> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition applies an “embedded” surface tension source term when
solving capillary hydrodynamics problems with phase function level set interface
tracking. Note that its counterpart for the base level set field is LS_CAPILLARY, and
this boundary condition is applied the same way to other level-set fields defined by the
Phase Function cards. It can be used with only subgrid integration. The surface
tension value used in this boundary condition is obtained from the Surface Tension
material parameter defined in the mat file. Note that each phase-function field requires
a separate PF_CAPILLARY boundary condition.

A description of the input parameters follows:

================ =================================================================
**PF_CAPILLARY** Name of the boundary condition.
**PF**           This string is used to indicated that this is a “boundary”
                 condition is applied at an internal phase boundary defined
                 by the zero contour of the level set function.
<integer>        An integer parameter that is used to specify to which phase
                 function field that the boundary condition is to be applied.
<float1>         Not currently used.
<float2>         Not currently used.
<float3>         Not currently used.
================ =================================================================

------------
**Examples**
------------

An example:
::

   BC = PF_CAPILLARY PF 1


-------------------------
**Technical Discussion**
-------------------------

Surface tension forces at a level set (phase function) representation of an interfacial
boundary are applied solely via this boundary condition. An additional divergence of
stress tensor term :math:`\Delta` ⋅ 
:math:`T_{cap}` is added to the fluid momentum equation. The form of this
tensor is

.. figure:: /figures/193_goma_physics.png
	:align: center
	:width: 90%

where *s* is the (isotropic) surface tension, *I* is the identity tensor, *n* is the vector normal
to the interface and :math:`\zeta_a` ( :math:`\phi` ) is the smoothed Dirac delta function. The surface tension
value used in this expression is obtained from the Surface Tension card found in the
material file.

The actual implementation in Goma integrates the divergence term by parts so the
expression that is added to the weak form of the momentum equation is:

.. figure:: /figures/194_goma_physics.png
	:align: center
	:width: 90%

This fact introduces the issue of integration error into the problem. As obvious above,
this source term involves the non-linear Dirac delta function factor. Conventional
numerical integration methods often do not offer adequate accuracy in evaluating this
integral, especially if if the interface width is a fraction of the average element size.
This has led to introduction the level-set-specific integration methods: subelement
integration and subgrid integration. In the latter case, more integration points are
clustered around the interface (in essence) to improve accuracy. The integer parameter
on the card should be set to zero to signify that the surface tension forces are distributed
in equal measure on both sides of the interfacial curve.

In the subelement integration case, however, an actual subelement mesh is place on
each of the interface-containing elements which is made to conform to the interface
curve. That is, the interface curve itself is covered by these subelement boundaries.
This allows the volume integral to be collapsed into a line integral and the line integral
evaluated along the subelement boundaries. This, however, introduces the problem of
identifying which side of the element the surface tension forces should actually be
applied to. Applying them to both simultaneously while either result in a cancellation
or a doubling of the surface tension effect. For these cases, the integer parameter on
this card is set to a -1 or a +1 to signify that the surface tension forces are applied to the
negative or positive side of the interface curve, respectively.



--------------
**References**
--------------

No References. 

..
	TODO - Lines 58 and 70 have pictures that need to be exhanged with equations.