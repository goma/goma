********************
**KIN_DISPLACEMENT**
********************

::

	BC = KIN_DISPLACEMENT SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition card is used as a distinguishing condition on the mesh motion
equations (viz. *mesh1, mesh2*, and *mesh3* under the *EQ* card). It forces the boundary of
the mesh defined by the side set to conform to a transient or steady material surface.
Unlike the *KINEMATIC* condition, which is designed for material surfaces between
two fluids, or the external material boundary of a fluid, this condition is applied to solid
materials to which the *TOTAL_ALE* mesh motion scheme is applied (see technical
discussion below and the *Mesh Motion* card). In two dimensions, this condition is
automatically applied to the normal component of the vector mesh equations, which is
rotated into normal-tangential form. In three dimensions, the application of this
boundary condition needs to be further directed with the *ROT* cards (see *ROTATION
specifications). The application of this condition should be compared with
KIN_DISPLACEMENT_PETROV* and *KIN_DISPLACEMENT_COLLOC*.

Definitions of the input parameters are as follows:

==================== ===============================================================
**KIN_DISPLACEMENT** Name of the boundary condition (<bc_name>).
**SS**               Type of boundary condition (<bc_type>), where SS denotes
                     side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side set in
                     EXODUS II) in the problem domain.
<integer>            Element block identification number for the region of TALE
                     solid mesh motion.
==================== ===============================================================

------------
**Examples**
------------

The following sample card:
::

     BC = KIN_DISPLACEMENT SS 7 12

leads to the application of the kinematic boundary condition (displacement form, see
below) to the boundary-normal component of the mesh-stress equation to the boundary
defined by side set 7. The element block ID number which shares this boundary with a
neighboring *TALE* or fluid *ARBITRARY* region is 12.

-------------------------
**Technical Discussion**
-------------------------

The functional form of the kinematic boundary condition is:


.. math::

    \underline n \cdot \left(\underline d_m - \underline d_m^0 \right) - \underline n \cdot \left(\underline d - \underline d^0 \right) = 0



Here **EQUATION** is the unit normal vector to the solid-fluid free surface, **EQUATION** 
is the mesh
displacement at the boundary **EQUATION**, is the mesh displacement from the base reference state
(which is automatically updated from the stress-free state coordinates and for
remeshes, etc. in *Goma* and need not be specified), **EQUATION** is the real solid displacement,
and **EQUATION** is the real solid displacement from the base reference state (or mesh). In stark
contrast with the *KINEMATIC* condition, which too is used to distinguish a material
fluid surface) this condition is written in Lagrangian displacement variables for *TALE*
mesh motion and is applied as a distinguishing condition on the mesh between a fluid
and *TALE* solid region. In essence, it maintains a real solid displacement field such that
no real-solid mass penetrates the boundary described by this condition.



--------------
**References**
--------------

SAND2000-0807: TALE: An Arbitrary Lagrangian-Eulerian Approach to Fluid-
Structure Interaction Problems, P. R. Schunk, May 2000

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

..
	 TODO - The picture in line 61 needs to be exchanged with the equation. In lines 65-75, where it says "**EQUATION**" there is supposed to be something from the equation that needs to be written. 