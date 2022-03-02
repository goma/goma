***************
**VELO_NORMAL**
***************

::

	BC = VELO_NORMAL SS <bc_id> <float> [integer]

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition allows the user to set the outward velocity component normal
to a surface.

Definitions of the input parameters are as follows:

================ ================================================================
**VELO_NORMAL**  Boundary condition designation.
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<float>          :math:`V_n`, value of the normal velocity component. Note that
                 this velocity component is relative to the motion of the
                 underlying mesh.
[integer]        *blk_id*, an optional parameter that is the element block
                 number in conjugate problems that identifies the
                 material region where the *VELO_NORMAL* condition
                 will be applied (usually the liquid element block in
                 solid/liquid conjugate problems). For external
                 boundaries, this optional parameter can be set to unity to
                 force the condition to be kept at a corner between two
                 side sets (2D only). This is handy for corner conditions.
                 Please see GTM-004.0 for details.
================ ================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_NORMAL SS 10   0.0

This boundary condition will enforce an impenetrability constraint over side set 10 as it
excludes normal velocity of the fluid relative to the mesh. This is by far the most
common context for this boundary condition.

-------------------------
**Technical Discussion**
-------------------------

* The actual weighted residual equation that is applied to a node, *j*, on the surface in
  question is as follows:

.. figure:: /figures/076_goma_physics.png
	:align: center
	:width: 90%

where :math:`\phi_j` is the finite element trial function, *n* the outward-pointing normal to the
surface, *v* the fluid velocity, :math:`v_s` the velocity of the underlying mesh, and :math:`v_n`
is the
normal velocity set by :math:`V_n` (the input value).

* This constraint is a rotated strongly integrated equation so that it will replace one
  of the rotated components of the fluid momentum equation. This component
  should generally always be the normal rotated component. In two dimensions, this
  replacement is automatic. In three dimensions, this replacement must be specified
  by a *ROT* condition.

* This card applies the identical constraint that is applied by the *KINEMATIC*
  boundary condition. The only difference is that this card replaces the normal
  component of the rotated *fluid momentum* equation, while the latter card replaces
  the normal component of the rotated (*pseudo-solid*) *mesh momentum* equation.

* In conjugate liquid/solid problems, the *VELO_NORMAL* condition is often used to
  enforce the impenetrability condition of the liquid/solid interface. The optional
  *blk_id* parameter can be used to insure that the *VELO_NORMAL* condition is
  correctly applied to the liquid side of the interface. *blk_id* should be set equal to the
  element block ID of the liquid in this case. This also applies to the *KINEMATIC*
  and *KINEMATIC_PETROV* boundary conditions.



--------------
**References**
--------------

GT-001.4: GOMA and SEAMS tutorial for new users, February 18, 2002, P. R. Schunk
and D. A. Labreche

GTM-004.1: Corners and Outflow Boundary Conditions in Goma, April 24, 2001, P. R.
Schunk