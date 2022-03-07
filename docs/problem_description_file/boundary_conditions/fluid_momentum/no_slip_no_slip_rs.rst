**********************
**NO_SLIP/NO_SLIP_RS**
**********************

::

	BC = {NO_SLIP | NO_SLIP_RS} SS <bc_id> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ VECTOR MOMENTUM)**

This card invokes a special boundary condition that applies a no-slip condition to the
fluid velocity at an interface between a liquid phase and a solid phase so that the fluid
velocity and solid velocity will be in concert. The solid phase must be treated as a
Lagrangian solid and may be in a convected frame of reference. The fluid velocity is
equal to the velocity of the stress-free state mapped into the deformed state (for steadystate
problems).

In general, a *SOLID_FLUID* boundary condition must also be applied to the same
boundary so that the force balance between liquid and solid is enforced. Note that a
*FLUID_SOLID* boundary condition will have no effect since the strongly enforced
*NO_SLIP/NO_SLIP_RS* on the fluid momentum equation will clobber it.

All elements on both sides of the interface must have the same element type, i.e., the
same order of interpolation and basis functions, e.g., Q1 or Q2.

Definitions of the input parameters are as follows:

+----------------------------+-------------------------------------------------------------+
|**{NO_SLIP | NO_SLIP_RS}**  | Boundary condition name applied in the following            |
|                            | formulations:                                               |
|                            |                                                             |
|                            |   * **NO_SLIP** - this condition applies when the solid     |
|                            |     phase is a purely *LAGRANGIAN* solid.                   |
|                            |   * **NO_SLIP_RS** - this condition should be used instead  |
|                            |     when the displacements in the solid phase are           |
|                            |     determined via a *TALE* formulation.                    |
+----------------------------+-------------------------------------------------------------+
|**SS**                      | Type of boundary condition (<bc_type>), where **SS**        |
|                            | denotes side set in the EXODUS II database.                 |
+----------------------------+-------------------------------------------------------------+
|<bc_id>                     | The boundary flag identifier, an integer associated with    |
|                            | <bc_type> that identifies the boundary location (side set   |
|                            | in EXODUS II) in the problem domain. This side set          |
|                            | should be the intersection of liquid and solid element      |
|                            | blocks and be defined so that it is present in both         |
|                            | element blocks.                                             |
+----------------------------+-------------------------------------------------------------+
|<integer1>                  | the element block ID number of the solid phase material.    |
+----------------------------+-------------------------------------------------------------+
|<integer2>                  | the element block ID number of the liquid phase material.   |
+----------------------------+-------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:
::

     BC= NO_SLIP SS 10   2 1

This card will enforce continuity of velocity between the solid phase in element block 2
with the fluid phase in element block 1. Side set 10 should be in common with both
element blocks.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition is a vector condition meaning that all three components
  of the fluid momentum equation are affected by use of a single boundary
  condition. The actual constraints that are imposed at node *j* are:

.. figure:: /figures/075_goma_physics.png
	:align: center
	:width: 90%


where :math:`\phi_j` is the finite element trial function, :math:`v_f` is the fluid velocity, and :math:`v_s` is the
solid phase velocity. These three constraints are strongly enforced so they
replace completely the x, y, and z fluid momentum components. The boundary
condition is not rotated since all three components of the momentum equation
are supplanted.

* As mentioned above this boundary condition is used most often in conjunction
  with the *SOLID_FLUID* boundary condition which equates stresses across fluid/
  solid interfaces. As described in the section discussing this card, this latter card
  imposes these forces by using the residuals of the fluid momentum equation as
  surrogates for the fluid phase forces. These forces however are imposed on the
  solid equations prior to imposition of the *NO_SLIP* boundary condition.

* As noted above, for this boundary condition to function properly it is necessary
  that the side set between the fluid and solid element block be present in both
  element blocks. To explain this it is necessary to recognize that side sets are
  defined as a set of faces attached to specific elements. This is in contrast to node
  sets which are simply a list of node numbers. Therefore, in the case of a side set
  that lies at the interface of two element blocks, it is possible for a given face in that 
  side set to appear twice, once attached to the element in the first element block and
  a second time attached to the adjoining element in the second element block. This
  is the condition that is required for the proper execution of this boundary
  condition. Fortunately, this is the default of most meshing tools that interface with
  *Goma*.

* It is also important to reiterate that another necessary condition for the proper
  function of this boundary condition is that the interpolation order of the pseudosolid
  mesh unknowns *and the fluid velocity unknowns* in the *ALE* fluid phase block
  be identical to the interpolation order of the solid displacement unknowns in the
  *LAGRANGIAN* or *TALE* adjoining solid phase block. This usually means that the
  element type must be the same in both phases. In two-dimensions this generally is
  not a problem, but in three dimensions it can impose a considerable hardship on
  the analyst.



--------------
**References**
--------------

No References.