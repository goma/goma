*************
**MOVING_CA**
*************

::

	BC = MOVING_CA NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MOMENTUM)**

The intent of this boundary condition is to apply a contact angle at wetting in a twodimensional
flow that is a function of the rate of advance or recession of the contact
line over the substrate. It is experimental, untested, and unsupported; use it at your own
risk.

There are ten values that must be specified in the <float_list>; definitions of the input
parameters are as follows:

============== =================================================================
**MOVING_CA**  Name of the boundary condition.
**NS**         Type of boundary condition (<bc_type>), where **NS** denotes
               node set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (node set in
               EXODUS II) in the problem domain. A *KINEMATIC* free
               surface must terminate at this single-node node set.
<float1>       :math:`\theta_{stc}`, the static contact angle, in degrees.
<float2>       :math:`n_x`, the x-component of solid surface normal vector.
<float3>       :math:`n_y`, the y-component of solid surface normal vector.
<float4>       :math:`n_z`, the z-component of solid surface normal vector.
<float5>       :math:`\theta_{adv}`, the advancing contact angle, in degrees.
<float6>       :math:`\theta_{rec}`, the receding contact angle, in degrees.
<float7>       :math:`\alpha`, the scale-factor, see below.
<float8>       :math:`v_{wx}`, the x-component of wall velocity.
<float9>       :math:`v_{wy}`, the y-component of wall velocity.
<float10>      :math:`v_{wz}`, the z-component of wall velocity.
============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

    BC = MOVING_CA NS 100 90.0 0. 1. 0. 135.0 45.0 1.0 -1. 0. 0.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition applies a point collocated constraint on the angle between
  the solid surface vector and the free-surface normal of the form:

.. figure:: /figures/111_goma_physics.png
	:align: center
	:width: 90%

where *n* is the solid surface vector specified on the card and :math:`n_{fs}` is the free-surface
normal computed automatically by *Goma*. The contact angle is variable depending
upon the relative velocity of the mesh speed, :math:`\dot{x}` , and the substrate speed, :math:`v_w`
specified on the card float_list:

.. figure:: /figures/112_goma_physics.png
	:align: center
	:width: 90%

* This constraint on the moving contact angle replaces a rotated component of the
  momentum equation. In effect a wetting force is applied at the contact line whose
  magnitude depends on the discrepancy between actual contact angle and that
  computed by the above expressions. Note that other contact angle constraints are
  applied to rotated components of the mesh equation. A real question exists whether
  such a formulation is consistent with a *KINEMATIC* boundary condition also
  applied to this node.

* Not also that since this boundary condition is applied to the momentum equation,
  care must be taken to relax any Dirichlet on the substrate velocity. Otherwise, this
  latter constraint will override this constraint.

* Users are again cautioned that this boundary condition is untested and potentially
  inconsistent. It may not work.




.. TODO - Lines 59 and 68 have photos that needs to be replaced with the real equation.