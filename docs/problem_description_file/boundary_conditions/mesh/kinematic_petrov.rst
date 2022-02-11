********************
**KINEMATIC_PETROV**
********************

::

	BC = KINEMATIC_PETROV SS <bc_id> <float1> [integer]

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition on the mesh motion
equations (viz. *mesh1, mesh2*, and *mesh3* under the *EQ* card). It enforces the boundary
of the mesh defined by the side set to conform to a transient or steady material surface,
with an optional, pre-specified mass loss/gain rate. In two dimensions, this condition is
automatically applied to the normal component of the vector mesh equations, which is
rotated into normal-tangential form. In three dimensions, the application of this
boundary condition needs to be further directed with the *ROT* cards (see *ROTATION
Specifications*). Please consult the Technical Discussion for important inofrmation.

Definitions of the input parameters are as follows:

==================== ==============================================================
**KINEMATIC_PETROV** Name of the boundary condition (<bc_name>).
**SS**               Type of boundary condition (<bc_type>), where **SS** denotes
                     side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side set in
                     EXODUS II) in the problem domain.
<float1>             Mass-loss (positive) or mass-gain (negative) velocity at the
                     free boundary.
[integer]            Optional integer value indicating the element block id from
                     which to apply the boundary condition.
==================== ==============================================================

------------
**Examples**
------------

The following sample card
::

     BC = KINEMATIC_PETROV SS 7 0.0

leads to the application of the kinematic boundary condition to the boundary-normal
component of the mesh-stress equation to the boundary defined by side set 7.

-------------------------
**Technical Discussion**
-------------------------

**Important note**: This condition is actually the same as the *KINEMATIC* condition but
is applied with different numerics for special cases. Specifically, rather than treated in a
Galerkin fashion with a weighting function equal to the interpolation function for
velocity, the residual of the equation is formed as weighted by the directional derivative
of the basis functions along the free surface. Specifically,

.. math::

   \int \left\{ \underline{n} \cdot \left( \underline{v} - \underline{v}_s \right) - \dot{m} \right\} \phi^i \, \mathrm{d}A = R^i = 0

where the nodal basis function :math:`\phi^i` is replaced by :math:`\frac{\partial}{\partial s} \phi^i` in the residual equation. Compare
this to the *KINEMATIC* boundary condition description.




.. 
	TODO - The picture in line 61 needs to be replaced with the equation and in line 65 where it says, "**EQUATION**" it needs to be written out. 
