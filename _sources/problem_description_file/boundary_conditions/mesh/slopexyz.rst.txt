************
**SLOPEXYZ**
************

::

	BC = {SLOPEX | SLOPEY | SLOPEZ} SS <bc_id> <floatlist>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MESH)**

This boundary condition card applies a slope at the boundary of a *LAGRANGIAN,
TALE*, or *ARBITRARY* solid (see *Mesh Motion* card) such that the normal vector to the
surface is colinear with the vector specified as input, viz 
:math:`\underline{n} \cdot \underline{n}_{\mathrm{spec}} = 0`. Here :math:`\underline{n}_{\mathrm{spec}}`
is the vector specified component-wise via the three <floatlist> parameters on the input
card. Definitions of the input parameters are as follows:

=============================== ==============================================================
**{SLOPEX | SLOPEY | SLOPEZ}**  Name of the boundary condition (<bc_name>).
**SS**                          Type of boundary condition (<bc_type>), where **SS** denotes
                                side set in the EXODUS II database
<bc_id>                         The boundary flag identifier, an integer associated with
                                <bc_type> that identifies the boundary location (side set in
                                EXODUS II) in the problem domain.
<float1>                        X-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
<float2>                        Y-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
<float3>                        Z-component of the slope vector :math:`\underline{n}_{\mathrm{spec}}`
=============================== ==============================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = SLOPEX SS 10 1.0 1.0 0.0

This card invokes a boundary condition on the normal component of the mesh residual
momentum equations such that the outward facing surface normal vector along side set
10 is colinear with the vector [1.0, 1.0, 0.0]. This condition is applied to the x-component
of the mesh residual equations.

-------------------------
**Technical Discussion**
-------------------------

See discussion for *BC* card *SLOPE*. The only difference in these conditions and the
*SLOPE* conditions, is that the latter invokes rotation of the vector mesh residual
equations on the boundary.



