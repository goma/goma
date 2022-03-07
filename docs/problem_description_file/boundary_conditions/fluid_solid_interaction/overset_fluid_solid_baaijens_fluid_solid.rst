********************************************
**OVERSET_FLUID_SOLID/BAAIJENS_FLUID_SOLID**
********************************************

::

	BC = BAAIJENS_FLUID_SOLID PF <pf_id> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(EMBEDDED_SURF/R_MOMENTUM1)**

This boundary condition is used to apply a traction to the fluid that comes from a solid
while using *Goma’s* overset grid capability. The condition is used when the complete
fluid-structure interaction problem is being solved, viz. stresses between fluid and solid
are both accommodated as is the dynamics of the structure and fluid. The condition is
applied to the fluid phase along a zero-level-set contour, hence the PF BC ID type. In
another mode of usage, Goma allows for a structure to be moved through a fluid under
prescribed kinematics, and in that case this condition is still applied as a solid traction
to the fluid. The value of that traction is dictated by the Lagrange multiplier kinematic
constraint (cf. LANGRANGE_NO_SLIP BC and LS_NO_SLIP BC). Note that the
condition is applied to a boundary in the fluid defined by a phase-field function (see
*phase1* equation type). . Two integer inputs together with a sideset ID integer are
required for this boundary condition:

======================= =================================================================
**BAAIJENS_FLUID_SOLD** Name of the boundary condition.
**PF**                  Type of boundary condition (<bc_type>), where **PF** denotes
                        a surface defined by a phase function (level-set).
<pf_id>                 The boundary flag identifier basically sets the number of the
                        phase field function to which this condition applies. For
                        now you must set this to 1, as this phase-field is hardwired
                        to handle the imprinted fluid solid boundary.
<integer1>              Element block ID of solid phase from the EXODUS II
                        database.
<integer2>              Element block ID of liquid phase from the EXODUS II
                        database.
======================= =================================================================

The peculiar name was derived from a paper by Frank Baaijens, from which
Goma’s formulation was generated. We are in the process of changing that
name to OVERSET_FLUID_SOLID.

------------
**Examples**
------------

Following is a sample card:
::

   BC = BAAIJENS_FLUID_SOLID PF 1 1 2

::

   BC = LS_NO_SLIP PF 1 1 2

This condition set applies a fluid traction condition to a surface defined by phase field
1, which is slaved to a side set that is set in the phase function slave surface capability.
(see Phase Function Initialization Method).

-------------------------
**Technical Discussion**
-------------------------

See discussion on LANGRANGE_NO_SLIP. This condition applies the fluid traction
boundary term on the fluid momentum equation.



--------------
**References**
--------------

GT-026.3