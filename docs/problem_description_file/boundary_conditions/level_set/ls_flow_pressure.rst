********************
**LS_FLOW_PRESSURE**
********************

::

	BC = LS_FLOW_PRESSURE LS <integer> <float1>

-----------------------
**Description / Usage**
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition applies a scalar pressure value as an “embedded” source term
on the fluid momentum equation at the zero level set contour. It can be used both when
subgrid or subelement integration is being used.

A description of the input parameters follows:

===================== ===============================================================
**LS_FLOW_PRESSURE**  Name of the boundary condition.
**LS**                This string is used to indicated that this is a “boundary”
                      condition is applied at an internal phase boundary defined
                      by the zero contour of the level set function.
<integer>             An integer parameter than is permitted to take one of three
                      values -1, 0, or 1. Depending upon the choice of this
                      parameter the pressure value is applied to the negative
                      phase, both phase, or the positive phase, respectively.
                      Details are given below.
<float1>              *P*, The constant value of pressure to be applied at the zero
                      level set contour.
===================== ===============================================================

------------
**Examples**
------------

An example:
::

   BC = LS_FLOW_PRESSURE LS 0 1013250.0

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is somewhat analogous to the FLOW_PRESSURE boundary
condition used quite often in ALE problems. It applies a scalar pressure at the
interfacial curve as an embedded boundary condtion. It can be used in by subgrid and
subelement methods. In the case of the former, a distributed volume integral of the
form:

.. figure:: /figures/201_goma_physics.png
	:align: center
	:width: 90%

where :math:`\vec{n}_{fs}` is the normal to the level set contour 
:math:`\delta_\alpha` (:math:`\phi`) and is the familiar smoothed
Dirac delta function with width parameter :math:`\alpha`. When subelement integration is used this
width parameter goes to zero and the volume integral becomes a surface integral along
the zero level set contour (Note: as of Oct 2005 subelement integration is not
supported for three dimensional problems).

When using this boundary condition concurrent with subgrid integration, the integer
parameter that appears on the card should be consistently set to zero. This ensures the
volume source will be applied symmetrically. However, when using subelement
integration this integer parameter must be entire a +1 or a -1 so that the pressure force
will be applied to only on side of the interface and not both which would result in
cancellation. This is much the same as was seen for the LS_CAPILLARY boundary
condition and the reader is referred to that card for a more detailed discussion.



--------------
**References**
--------------

No References. 

..
	TODO - Line 54 has a picture that needs to be exhanged with equation.
