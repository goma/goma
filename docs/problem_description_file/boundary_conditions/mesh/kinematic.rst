*************
**KINEMATIC**
*************

::

	BC = KINEMATIC SS <bc_id> <float1> [integer]

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
boundary condition needs to be further directed with the *ROT* cards (see *Rotation
Specifications*). The application of this condition should be compared with
*KINEMATIC_PETROV* and *KINEMATIC_COLLOC*.

Definitions of the input parameters are as follows:

=============== ================================================================
**KINEMATIC**   Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS** denotes
                side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set in
                EXODUS II) in the problem domain.
<float1>        Mass-loss (positive) or mass-gain (negative) velocity at the
                free boundary.
[integer]       Optional integer value indicating the element block id from
                which to apply the boundary condition.
=============== ================================================================

------------
**Examples**
------------

The following sample card
::

     BC = KINEMATIC SS 7 0.0

leads to the application of the kinematic boundary condition to the boundary-normal
component of the mesh-stress equation on the boundary defined by side set 7.

-------------------------
**Technical Discussion**
-------------------------

The functional form of the kinematic boundary condition is:

.. math::

   \underline{n} \cdot \left( \underline{v} - \underline{v}_s \right) = \dot{m}

Here :math:`\underline{n}` is the unit normal vector to the free surface, :math:`\underline{v}` is the velocity of the fluid, :math:`\underline{v}_s` is
the velocity of the surface (or mesh), and :math:`\dot{m}` is the mass loss/gain rate. In two
dimensions this equation is applied to the normal component of the vector mesh
position equation, and hence is considered as a distinguishing condition on the location
of the mesh relative to the fluid domain.


--------
**FAQs**
--------

See the FAQ pertaining to “Continuation Strategies for Free Surface Flows” on the
*DISTNG* boundary condition card.

--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

.. 
	TODO - The image in line 58 needs to be changed to the equation. In lines 62-66 there are symbols for the equation above that need to be checked.