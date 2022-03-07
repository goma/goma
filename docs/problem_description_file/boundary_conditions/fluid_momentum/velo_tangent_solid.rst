**********************
**VELO_TANGENT_SOLID**
**********************

::

	BC = VELO_TANGENT_SOLID SS <bc_id> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition sets the tangential fluid velocity component at a fluid/solid
interface to the tangential velocity component of the solid material. The latter includes
any motion of the stress-free state. This boundary condition is applicable only to twodimensional
problems and is normally used in conjunction with the Total Arbitrary
Lagrangian/Eulerian algorithm in Goma (See GT-005.3).

Definitions of the input parameters are as follows:

====================== ========================================================
**VELO_TANGENT_SOLID** The name of the boundary condition.
**SS**                 Type of boundary condition (<bc_type>), where **SS**
                       denotes side set in the EXODUS II database.
<bc_id>                The boundary flag identifier, an integer associated with
                       <bc_type> that identifies the boundary location (side set
                       in EXODUS II) in the problem domain.
<integer1>             The element block id defining the solid phase adjacent
                       to <bc_id>.
<integer2>             The element block id defining the liquid phase adjacent
                       to <bc_id>.
====================== ========================================================

------------
**Examples**
------------

The following is an example of this card
::

     BC = VELO_TANGENT_SOLID SS 10   2 1

In this case, sideset 10 is an internal sideset between two separate materials, the solid
material in element block 2 and the liquid material in element block 1.

-------------------------
**Technical Discussion**
-------------------------

The boundary condition being applied is the strong integrated condition:

.. figure:: /figures/092_goma_physics.png
	:align: center
	:width: 90%

where :math:`v_m` is the fluid velocity, :math:`v_{sfs}` is the velocity of the solid material stress-free-state
(usually solid-body translation, or rotation..see *Advected Langragian Velocity* card)
including the motion of the deformed coordinates, and *t* is the vector tangent to the side
set. :math:`F_m` is the deformation gradient tensor and the time derivative term is the motion of
the deformed state tangential to the surface in question.

This condition is advocated for use with the TALE algorithm (see GT-005.3).




.. TODO - Line 54 contains a photo that needs to be exchanged for the equation.