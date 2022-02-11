************
**LS_INLET**
************

::

	BC = LS_INLET SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(PCC/LEVEL SET)**

This boundary condition is used to set the values of the level set function on a sideset.
Most of this is done on an inlet or outlet boundary to elminate the potential for
oscillations in the level set field at those points from introducing spurious interfacial
(zero contours) .

A description of the input parameters follows:

============ ===================================================================
**LS_INLET** Name of the boundary condition.
**SS**       This string indicates what this boundary is applied to.
<bc_id>      An integer parameter than is permitted to take one of three
             values -1, 0, or 1. Depending upon the choice of this
             parameter the surface tension forces are applied to the
             negative phase, both phase, or the positive phase,
             respectively. Details are given below.
============ ===================================================================

------------
**Examples**
------------

An example:
::

   BC = LS_INLET SS 10

-------------------------
**Technical Discussion**
-------------------------

Surface tension forces at a level set representation of an interfacial boundary are
applied solely via this boundary condition. An additional divergence of stress tensor
term :math:`\Delta` ⋅ :math:`T_{cap}` is added to the fluid momentum equation. Following Jacqmin the form of this tensor is

.. figure:: /figures/202_goma_physics.png
	:align: center
	:width: 90%

where s is the (isotropic) surface tension, I is the identity tensor, n is the vector normal to the interface and :math:`\delta_\alpha` (:math:`\phi`) is the smoothed Dirac delta function. The surface tension
value used in this expression is obtained from the Surface Tension card found in the
material file.

The actual implementation in Goma integrates the divergence term by parts so the
expression that is added to the weak form of the momentum equation is:

.. figure:: /figures/203_goma_physics.png
	:align: center
	:width: 90%

This fact introduces the issue of integration error into the problem. As obvious above,
this source term involves the non-linear Dirac delta function factor. Conventional
numerical integration methods often do not offer adequate accuracy in evaluating this
integral, especially if if the interface width is a fraction of the average element size.
This has led to introduction the level-set-specific integration methods: subelement
integration and subgrid integration. In the latter case, more integration points are
clustered around the interface (in essence) to improve accuracy. The integer parameter
on the card should be set to zero to signify that the surface tension forces are 
distributed in equal measure on both sides of the interfacial curve.

In the subelement integration case, however, an actual subelement mesh is place on
each of the interface-containing elements which is made to conform to the interface
curve. That is, the interface curve itself is covered by these subelement boundaries.
This allows the volume integral to be collapsed into a line integral and the line integral
evaluated along the subelement boundaries. This, however, introduces the problem of
identifying which side of the element the surface tension forces should actually be
applied to. Applying them to both simultaneously while either result in a cancellation
or a doubling of the surface tension effect. For these cases, the integer parameter on
this card is set to a -1 or a +1 to signify that the surface tension forces are applied 
to the negative or positive side of the interface curve, respectively.



--------------
**References**
--------------

No References. 

..
	TODO - Lines 50 and 61 have pictures that needs to be exhanged with equations.
