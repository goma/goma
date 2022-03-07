****************
**LS_CAPILLARY**
****************

::

	BC = LS_CAPILLARY LS <integer>

-----------------------
**Description / Usage**
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition applies an “embedded” surface tension source term when
solving capillary hydrodynamics problems with level set interface tracking. It can be
used both when subgrid or subelement integration is being used. The surface tension
value used in this boundary condition is obtained from the Surface Tension material
parameter defined in the mat file.

A description of the input parameters follows:

================= =======================================================
**LS_CAPILLARY**  Name of the boundary condition.
**LS**            This string is used to indicated that this is a “boundary”
                  condition is applied at an internal phase boundary defined
                  by the zero contour of the level set function.
<integer>         An integer parameter than is permitted to take one of three
                  values -1, 0, or 1. Depending upon the choice of this
                  parameter the surface tension forces are applied to the
                  negative phase, both phase, or the positive phase,
                  respectively. Details are given below.
================= =======================================================

------------
**Examples**
------------

An example:
::

   BC = LS_CAPILLARY LS 0

-------------------------
**Technical Discussion**
-------------------------

*First, a warning: If subelement integration is off, make sure there is a nonzero levelset
length scale or no surface forces term will be applied.*

Surface tension forces at a level set representation of an interfacial boundary are
applied solely via this boundary condition. An additional divergence of stress tensor
term :math:`\Delta` ⋅ :math:`T_{cap}` is added to the fluid momentum equation. Following Jacqmin the
form of this tensor is

.. figure:: /figures/199_goma_physics.png
	:align: center
	:width: 90%

where s is the (isotropic) surface tension, I is the identity tensor, *n* is the vector normal
to the interface and :math:`\zeta_a` ( :math:`\phi` ) is the smoothed Dirac delta function. The surface tension
value used in this expression is obtained from the Surface Tension card found in the
material file.

The actual implementation in Goma integrates the divergence term by parts so the
expression that is added to the weak form of the momentum equation is:

.. figure:: /figures/200_goma_physics.png
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
	TODO - Lines 56 and 68 have pictures that need to be exhanged with equations.
