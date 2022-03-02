************
**LS_YFLUX**
************

::

	BC = LS_YFLUX LS <integer> <integer1> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

**(EMB/ENERGY)**

This boundary condition applies a scalar mass flux value as an “embedded” source
term on a species conservation equation at the zero level set contour. It can be used
both when subgrid or subelement integration is being used.

A description of the input parameters follows:

============= ==============================================================
**LS_YFLUX**  Name of the boundary condition.
**LS**        This string is used to indicated that this is a “boundary”
              condition is applied at an internal phase boundary defined
              by the zero contour of the level set function.
<integer>     An integer parameter than is permitted to take one of three
              values -1, 0, or 1. Depending upon the choice of this
              parameter the mass flux value is applied to the negative
              phase, both phase, or the positive phase, respectively.
              Details are given below.
<integer1>    *w*, This the species equation index to which this boundary
              condition is applied.
<float1>      :math:`h_c`, a constant value for the mass transfer 
              coefficient at the interface. 
<float2>      :math:`Y_c`, the “bulk” concentration of species used in 
              conjunction with the mass transfer coefficient to compute
              the mass flux.
============= ==============================================================

------------
**Examples**
------------

An example:
::

   BC = LS_YFLUX LS 0 0 1.e-2 0.75

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is somewhat analogous to the YFLUX boundary condition
used quite often in non-level set problems to apply a scalar species flux at a at boundary
defined by a side set. It applies a scalar mass transfer flux at the interfacial curve as an
embedded boundary condtion. It can be used in by subgrid and subelement methods.
In the case of the former, a distributed volume integral of the form:

.. figure:: /figures/207_goma_physics.png
	:align: center
	:width: 90%

where :math:`\delta_\alpha` (:math:`\phi`) is the familiar smoothed Dirac delta function with width parameter :math:`\alpha` and the mass flux , J, is given by the typical relation:

.. figure:: /figures/208_goma_physics.png
	:align: center
	:width: 90%

When subelement integration is used this width parameter goes to zero and the volume
integral becomes a surface integral along the zero level set contour (Note: as of Oct
2005 subelement integration is not supported for three dimensional problems).

When using this boundary condition concurrent with subgrid integration, the integer
parameter that appears on the card should be consistently set to zero. This ensures the
volume source will be applied symmetrically. However, when using subelement
integration this integer parameter must be entire a +1 or a -1 so that the mass flux will
be applied only on side of the interface and not both which would result in cancellation.
This is much the same as was seen for the LS_CAPILLARY boundary condition and
the reader is referred to that card for a more detailed discussion.



--------------
**References**
--------------

No References. 

.. TODO -Lines 59 and 64 have pictures that need to be swapped with the correct equations.