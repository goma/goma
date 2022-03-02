*********
**SLOPE**
*********

::

	BC = SLOPE SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card applies a slope at the boundary of a *LAGRANGIAN,
TALE*, or *ARBITRARY* solid (see *Mesh Motion* card) such that the normal vector to the
surface is colinear with the vector specified as input, viz :math:`\underline{n} \cdot \underline{n}_{\mathrm{spec}} = 0` . Here :math:`\underline{n}_{\mathrm{spec}}`
the vector specified component-wise via the three <float> parameters on the input card.
Definitions of the input parameters are as follows:

========== ===================================================================
**SLOPE**  Name of the boundary condition (<bc_name>).
**SS**     Type of boundary condition (<bc_type>), where **SS** denotes
           side set in the EXODUS II database
<bc_id>    The boundary flag identifier, an integer associated with
           <bc_type> that identifies the boundary location (side set in
           EXODUS II) in the problem domain.
<float1>   X-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
<float2>   Y-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
<float3>   Z-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
========== ===================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = SLOPE SS 10 1.0 1.0 0.0

This card invokes a boundary condition on the normal component of the mesh residual
momentum equations such that the outward facing surface normal vector along side set
10 is colinear with the vector [1.0, 1.0, 0.0].

-------------------------
**Technical Discussion**
-------------------------

This condition, although not often used, allows for a planar boundary condition (cf.
*PLANE, PLANEX*, etc.) to be specified in terms of a slope, rather than a specific
equation. Clearly, at some point along the surface (most likely at the ends), the
geometry has to be pinned with some other boundary condition (cf. *DX, DY, DZ*) so as
to make the equation unique. This condition has the following mathematical form:

.. math::
   \underline{n} \cdot \underline{n}_{\mathrm{spec}} = 0

and is applied in place of the normal component of the mesh motion equations, i.e., it is
a rotated type boundary condition. If used in three dimensions, it will require a rotation
description with the *ROT* cards.



