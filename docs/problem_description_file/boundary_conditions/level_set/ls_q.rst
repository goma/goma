********
**LS_Q**
********

::

	BC = LS_Q LS <integer> <float1>

-----------------------
**Description / Usage**
-----------------------

**(EMB/ENERGY)**

This boundary condition applies a scalar heat flux value as an “embedded” source term
on the heat conservation equation at the zero level set contour. It can be used both
when subgrid or subelement integration is being used.

A description of the input parameters follows:

=========== ===============================================================
**LS_Q**    Name of the boundary condition.
**LS**      This string is used to indicated that this is a “boundary”
            condition is applied at an internal phase boundary defined
            by the zero contour of the level set function.
<integer>   An integer parameter than is permitted to take one of three
            values -1, 0, or 1. Depending upon the choice of this
            parameter the heat flux value is applied to the negative
            phase, both phase, or the positive phase, respectively.
            Details are given below.
<float1>    *q*, The constant value of heat flux to be applied at the zero
            level set contour.
=========== ===============================================================

------------
**Examples**
------------

An example:
::

   BC = LS_Q LS 0 -1.1e-3

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is somewhat analogous to the QSIDE boundary condition
used quite often in non-level set problems. It applies a scalar heat flux at the interfacial
curve as an embedded boundary condtion. It can be used in by subgrid and subelement
methods. In the case of the former, a distributed volume integral of the form:

.. figure:: /figures/206_goma_physics.png
	:align: center
	:width: 90%

where :math:`\Delta` ⋅ :math:`T_{cap}` is the familiar smoothed Dirac delta function with width parameter :math:`\alpha`.
When subelement integration is used this width parameter goes to zero and the volume
integral becomes a surface integral along the zero level set contour (Note: as of Oct
2005 subelement integration is not supported for three dimensional problems).

When using this boundary condition concurrent with subgrid integration, the integer
parameter that appears on the card should be consistently set to zero. This ensures the
volume source will be applied symmetrically. However, when using subelement
integration this integer parameter must be entire a +1 or a -1 so that the heat flux will be
applied only on side of the interface and not both which would result in cancellation.
This is much the same as was seen for the LS_CAPILLARY boundary condition and
the reader is referred to that card for a more detailed discussion.



--------------
**References**
--------------

No References. 

..
	TODO - Line 53 has a picture that needs to be exhanged with an equation.
