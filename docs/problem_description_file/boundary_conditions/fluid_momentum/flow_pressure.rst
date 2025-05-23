*****************
**FLOW_PRESSURE**
*****************

::

	BC = FLOW_PRESSURE SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition card is used to set a constant value of pressure on a boundary.
Most often this condition is used to set an upstream or downstream pressure over a
fully-developed inflow/outflow boundary.

Definitions of the input parameters are as follows:

================== ===================================================
**FLOW_PRESSURE**  Boundary condition name.
**SS**             Type of boundary condition (<bc_type>), where **SS**
                   denotes side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (side set
                   in EXODUS II) in the problem domain.
<float>            :math:`P_{ex}`, the applied pressure. Positive values imply
                   compressive forces on the fluid, negative values imply
                   tensile forces.
================== ===================================================

------------
**Examples**
------------

The following sample input card will impose a constant compressive pressure force on
the boundary defined by sideset 23:
::

     BC = FLOW_PRESSURE SS 23   5.0

-------------------------
**Technical Discussion**
-------------------------

* The actual boundary condition that is applied to the fluid is given as follows:

.. figure:: /figures/095_goma_physics.png
	:align: center
	:width: 90%

where *n* is the outward normal vector to the boundary, *T* is the total fluid stress
tensor, and *P* is the applied pressure equal to <float1> above. From this the user
should be able to deduce the appropriate sign for his/her pressure value.

* This boundary condition is a weak integrated condition implying that it is added to
  all three components of the fluid momentum equation prior to rotation of equations
  or application of strongly enforced conditions or Dirichlet conditions.

* The astute user who is also well-versed in finite element formulations and
  terminology will recognize that this boundary condition is providing a value for
  the boundary condition term that appears after application of the divergence
  theorem to the weighted fluid momentum residual equations. Hence, imposing a
  value of zero for <float1> is exactly equivalent to saying nothing at all about the
  fluid velocity at a boundary.

* This boundary condition is found predominantly in two applications. First, setting
  the external pressure imposed on a free surface, and second, providing the driving
  force for flow by being imposed on an inflow or outflow fully-developed
  boundary. In this latter role, the usual procedure is to apply the
  *FLOW_PRESSURE* condition while strongly enforcing a zero condition on the
  velocity components transverse to the boundary. For boundaries parallel to one of
  the principle coordinate directions, Dirichlet conditions can be used to set these
  transverse components. For other inflow or outflow boundaries, it is suggested that
  the *VELO_TANGENT* and *VELO_TANGENT_3D* cards be employed instead.

* This boundary condition is very useful when working with non-Newtonian models
  where the inlet velocity field is apt to be complicated and hard to determine *a*
  *priori*. By imposing a pressure at the inflow with this card, the non-Newtonian
  inlet velocity profile will be determined implicitly. Augmenting conditions can
  then be used to couple the imposed pressure to the average flow rate over the
  boundary for an even more advanced capability.




.. TODO - Line 49 contains a photo that needs to be exchanged for the equation.