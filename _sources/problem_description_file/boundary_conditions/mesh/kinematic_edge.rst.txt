******************
**KINEMATIC_EDGE**
******************

::

	BC = KINEMATIC_EDGE <bc_id1> <bc_id2> <float1>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition on the mesh motion
equations (viz. *mesh1, mesh2*, and *mesh3* under the *EQ* card). It enforces the boundary
of the mesh defined by the side set to conform to a transient or steady material surface,
with an optional, pre-specified mass loss/gain rate. This condition is applied only in
three-dimensional problems along contact lines that define the intersection of a freesurface
and a geometrical solid, the intersection of which is partially characterized by
the binormal tangent as described below.

Definitions of the input parameters are as follows:

================== ===================================================================
**KINEMATIC_EDGE** Name of the boundary condition (<bc_name>).
**SS**             Type of boundary condition (<bc_type>), where **SS**
                   denotes side set in the EXODUS II database.
<bc_id1>           The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set
                   in EXODUS II) in the problem domain. This surface is
                   the “primary solid surface”
<bc_id2>           The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set
                   in EXODUS II) in the problem domain. This surface is
                   the “free surface”
<float1>           Mass-loss (positive) or mass-gain (negative) velocity at
                   the free boundary.
================== ===================================================================

------------
**Examples**
------------

::

     BC = KINEMATIC_SPECIES SS 10 2 0.0

In this example, the *KINEMATIC_EDGE* boundary condition is applied to the line
defined by the intersection of side sets 10 and 20. The normal vector used in
application of this condition is the one in the plane of side-set 10, viz. it is tangent to
the surface delineated by side set 10.

-------------------------
**Technical Discussion**
-------------------------

The functional form of the kinematic boundary condition is:

.. math::

   \underline{n}_{\mathrm{cl}} \cdot \left( \underline{v} - \underline{v}_s \right) = 0

Here :math:`\underline{n}_{\mathrm{cl}}` is the unit normal tangent vector to a line in space defined by two surfaces, in
the plane of the primary surface, viz. tangent to that surface. :math:`\underline{v}` is the velocity of the
fluid, :math:`\underline{v}_s` is the velocity of the surface (or mesh). This condition only makes sense in
three dimensions, and needs to be directed with *ROT* conditions for proper application.

.. figure:: /figures/042_goma_physics.png
	:align: center
	:width: 90%



--------------
**References**
--------------

GT-007.2: Tutorial on droplet on incline problem, July 30, 1999, T. A. Baer