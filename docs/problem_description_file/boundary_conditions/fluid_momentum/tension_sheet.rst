*****************
**TENSION_SHEET**
*****************

::

	BC = TENSION_SHEET SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/VECTOR MOMENTUM)**

This boundary condition card is used to apply a membrane force to the fluid
momentum equation in order to model a membrane-like structure, viz. one with no
bending stiffness but with significant tension much larger than the fluid viscous
stresses. This boundary condition is basically the same mathematically as the capillary
condition, with the tension here specified instead of a capillary surface tension. The
only difference is the way in which it is applied: it is applied as a strong integrated
condition instead of a weak form condition.

Definitions of the input parameters are as follows:

================= ==================================================================
**TENSION_SHEET** Name of the boundary condition.
**SS**            Type of boundary condition (<bc_type>), where **SS** denotes
                  side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set in
                  EXODUS II) in the problem domain.
<float1>          Τ, Tension of the membrane.
================= ==================================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = TENSION_SHEET SS 12

-------------------------
**Technical Discussion**
-------------------------

Usage notes:

* Can only be applied in two dimensions.

* One caveat that must be mentioned is that the formulation of these two boundary
  conditions is not general and therefore they should only be applied to web geometries
  that are predominantly horizontal. That is, the x component of the normal vector to the
  web sideset should at each point be less than or equal to the y component.

* To set the slope of the membrane at an endpoint, see the BC = SHEET_ENDSLOPE card.

* This boundary condition that can be used to model the interaction of fluid with a thin
  sheet under a constant tension load. The sideset to which it is applied must be fully
  “wetted” by a fluid. Note that this boundary condition arises as a simplification of the
  tensioned-web shell equations (cf. shell_tension and shell_curvature equations as
  described in GT-033.0 and GT-027.0) subject to two simplifying assumptions:

  1. The sheet supports no bending moments. That is, it isn’t very rigid.

  2. The tension in the sheet is significantly larger than the viscous stresses in the fluid.

Given these assumptions this boundary condition can be used to model tensioned web
applications without having to resort to the shell equations. It is a strongly integrated,
rotated boundary condition on the mesh equations. It can only be used in twodimensional
applications.



--------------
**References**
--------------

GT-033.0 and GT-027.0.
