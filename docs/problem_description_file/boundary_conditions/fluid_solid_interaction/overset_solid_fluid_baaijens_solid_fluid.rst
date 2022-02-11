********************************************
**OVERSET_SOLID_FLUID/BAAIJENS_SOLID_FLUID**
********************************************

::

	BC = OVERSET_SOLID_FLUID SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(CONTACT_SURF/ MESH)**

This boundary condition is used to apply a traction to a solid that comes from the fluid
while using *Goma’s* overset grid capability. The condition is used when the complete
fluid-structure interaction problem is being solved, viz. stresses between fluid and solid
are both accommodated as is the dynamics of the structure and fluid. The condition is
applied to the solid phase along a side set that defines the fluid/solid interface. Two
integer inputs together with a sideset ID integer are required for this boundary
condition:

======================== ==============================================================
**BAAIJENS_SOLID_FLUID** Name of the boundary condition.
**SS**                   Type of boundary condition (<bc_type>), where **SS** denotes
                         side set in the EXODUS II database.
<ss_id>                  The boundary flag identifier that sets the side set number.
<integer1>               Element block ID of solid phase from the EXODUS II
                         database.
<integer2>               Element block ID of liquid phase from the EXODUS II
                         database.
======================== ==============================================================

The peculiar name was derived from a paper by Frank Baaijens, from which
Goma’s formulation was generated. We are in the process of changing that
name to OVERSET_SOLID_FLUID.

------------
**Examples**
------------

Following is a sample card set:
::

   BC = BAAIJENS_SOLID_FLUID SS 1 2 1

::

   BC = BAAIJENS_FLUID_SOLID PF 1 2 1

::

   BC = LAGRANGE_NO_SLIP SS 1 2 1

Here, the BAAIJENS_SOLID_FLUID cared applies a boundary fluid traction to a
solid phase defined by side set 1. In this case the solid phase material ID is 2 and the
fluid phase 1.

-------------------------
**Technical Discussion**
-------------------------

See discussion on LAGRANGE_NO_SLIP. Basically, this condition results in a
boundary traction set by the Lagrange multiplier constraint to be applied to the solid
momentum equation (note the weak term that appears on that equation).



--------------
**References**
--------------

GT-026.3
