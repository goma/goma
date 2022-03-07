********************
**KINEMATIC_COLLOC**
********************

::

	BC = KINEMATIC_COLLOC SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition on the mesh motion
equations (viz. *mesh1, mesh2*, and *mesh3* under the *EQ* card). It enforces the boundary
of the mesh defined by the side set to conform to a transient or steady material surface,
with an optional, pre-specified mass loss/gain rate. In two dimensions this condition is
automatically applied to the normal component of the vector mesh equations, which is
rotated into normal-tangential form. In three dimensions the application of this
boundary condition needs to be further directed with the *ROT* cards (see *Rotation
Specifications*). Definitions of the input parameters are as follows:

===================== =============================================================
**KINEMATIC_COLLOC**  Name of the boundary condition (<bc_name>).
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain.
<float1>              Mass-loss (positive) or mass-gain (negative) velocity at
                      the free boundary.
===================== =============================================================

------------
**Examples**
------------

The following sample card
::

     BC = KINEMATIC_COLLOC SS 7 0.0

leads to the application of the kinematic boundary condition to the boundary-normal
component of the mesh-stress equation to the boundary defined by side set 7.

-------------------------
**Technical Discussion**
-------------------------

**Important note**: This condition is actually the same as the *KINEMATIC* condition but
is applied with different numerics for special cases. Specifically, rather than treated in a
Galerkin fashion, with a weighting function equal to the interpolation function for
velocity, the residual equation is formed at each node directly, in a collocated fashion,
without Galerkin integration. This method is better suited for high-capillary number
cases in which Galerkin’s method is often not the best approach.



