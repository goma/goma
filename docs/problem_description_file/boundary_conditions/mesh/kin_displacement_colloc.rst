***************************
**KIN_DISPLACEMENT_COLLOC**
***************************

::

	BC = KIN_DISPLACEMENT_COLLOC SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

The *KIN_DISPLACEMENT_COLLOC* boundary condition is exactly the same as
*KIN_DISPLACEMENT* except in the way in which it is applied numerically to a
problem. See *KIN_DISPLACEMENT* for a full discussion.

Definitions of the input parameters are as follows:

=========================== ===============================================================
**KIN_DISPLACEMENT_COLLOC** Name of the boundary condition (<bc_name>).
**SS**                      Type of boundary condition (<bc_type>), where SS denotes
                            side set in the EXODUS II database.
<bc_id>                     The boundary flag identifier, an integer associated with
                            <bc_type> that identifies the boundary location (side set in
                            EXODUS II) in the problem domain.
<integer>                   Element block identification number for the region of TALE
                            solid mesh motion.
=========================== ===============================================================

Sometimes this condition is a better alternative to *KIN_DISPLACEMENT* to stabilize
the surface and prevent wiggles. If the user wants to know more regarding numerical
issues and implementation, consult the description for the fluid-counterpart
*KINEMATIC_COLLOC* card.

------------
**Examples**
------------

The following sample card:
::

     BC = KIN_DISPLACEMENT_COLLOC SS 7 12

leads to the application of the kinematic boundary condition (displacement form, see
below) to the boundary-normal component of the mesh-stress equation to the boundary
defined by side set 7. The element block ID number which shares this boundary with a
neighboring *TALE* or fluid *ARBITRARY* region is 12.

-------------------------
**Technical Discussion**
-------------------------

See discussions on the *KINEMATIC_COLLOC* and *KIN_DISPLACEMENT* cards.



--------------
**References**
--------------

No References.