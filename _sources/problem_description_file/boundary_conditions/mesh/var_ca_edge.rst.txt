***************
**VAR_CA_EDGE**
***************

::

	BC = VAR_ CA_EDGE SS <bc_id1> <bc_id2> <float_list>

-----------------------
**Description / Usage**
-----------------------

.. warning::

  The angle in this BC is specified in degrees, not radians. This is 
  different than the 2D equivalent.

**(SIC-EDGE/ROTATED MESH)**

This card is used to set a variable contact angle on a dynamic three-dimensional contact
line. A local contact angle is determined based upon the local rate of advance/recession
of the contact line with respect to a web, and is always associated with the secondary
sideset. This card specifies the static contact angle, θ\ :sub:`s`, and a linear proportionality
constant c\ :sub:`T` between the local advance/recession rate and the cosine of the contact
angle. The speed of the moving web is specified by components of the web velocity.
The contact angle is imposed between the outward-pointing normal of the primary
sideset and the outward-pointing normal of the secondary sideset.

Definitions of the input parameters are as follows:

================ ==============================================================
**VAR_CA_EDGE**  Name of the boundary condition.
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id1>         The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain. This identifies
                 the *primary side set*; it should be a free surface.
<bc_id2>         The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain. This identifies
                 the *secondary side set*, which should be a “fixed”
                 geometric entity, e.g. PLANE or SPLINE. Taken
                 together, the primary and secondary sidesets define an
                 edge over which this boundary is applicable.
<float1>         θ\ :sub:`s`, parameter that is the static contact angle, in degrees.
                 This is the contact angle that the fluid approaches when
                 the relative motion of the contact line and substrate is
                 zero.
<float2>         c\ :sub:`T`, parameter that is the linear proportionality constant
                 between the local advance/recession rate and the cosine
                 of the contact angle; see details below in the Technical
                 Discussion.
<float3>         W\ :sub:`x`, x-component of the substrate velocity.
<float4>         W\ :sub:`y`, y-component of the substrate velocity.
<float5>         W\ :sub:`z`, z-component of the substrate velocity.
================ ==============================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VAR_CA_EDGE SS 60 20   135. 0.02   0. -1. 0.

This card sets a variable contact angle condition on the edge between side sets 60 and
20. The static contact angle is 135 degrees and the slope parameter is 0.02. The solid
substrate is moving at the fixed velocity (0, -1., 0.).

-------------------------
**Technical Discussion**
-------------------------

* A contact line that moves relative to its underlying solid substrate is referred to as
  a dynamic contact line. For a dynamic contact line associated with threedimensional
  flows, it is recognized that the dynamic contact angle must change
  from point to point along the curve because the local advance/recession rate of the
  contact line *with respect to the substrate changes*. Taking this variability into
  account is the function of this card.

  To understand the function of this card, we first define that the advance (or
  recession) rate of the contact line, u\ :sub:`wet`, as the normal component of the contact
  line velocity, x\ :sub:`cl` , relative to the substrate velocity, *W*:

   .. math::

    u_{wet} = n_{cl} \cdot \left(W - \dot x_{cl} \right) 



  where n\ :sub:`cl` is a unit vector normal to the contact
  line in the plane of the substrate as illustrated
  in the sketch at right. For an advancing contact 
  line u\ :sub:`wet` is negative and the converse. We can
  also define a local capillary number by nondimensionalizing
  the advance rate as follows,

.. figure:: /figures/064_goma_physics.png
	:align: center
	:width: 90%

.. math::

  Ca_L = \mu u_{wet} / \sigma



|

  where μ is the viscosity and σ the surface tension.

.. figure:: /figures/066_goma_physics.png
	:align: center
	:width: 90%

|

  We choose to define the contact angle as the angle
  between the outward normal to the free-surface and
  the substrate normal pointing away from the fluid
  phase as illustrate here. From direct observation of
  contact lines, we know that increasing the advance
  rate will decrease the contact angle towards zero.
  Conversely, a decrease in the advance rate or increase
  of recession rate will increase the contact angle
  towards 180. We capture the essence of this behavior
  via a simple linear relationship between the local capillary number and the cosine
  of the contact angle:

.. math::

  cos\ \theta = cos\ \theta_s - c_T Ca_L



|

  where θ\ :sub:`s` and c\ :sub:`T` are two input parameters. The function of this card is to apply 
  this
  model for contact angle on the contact line curve.

* This model has many restrictions. It is really only valid for very very small | Ca\ :sub:`L` |
  and also does not predict that the contact angle asymptotically approaches 0 or 180
  for | Ca\ :sub:`L` | very large. Instead, it is algorithmically restricted to returning 0 or 180 if
  the above linear relation would predict an angle outside of these bounds.

* Unlike the *CA_EDGE* boundary condition, the *VAR_CA_EDGE* condition is
  applied as a strong integrated constraint. The equation associated with each node
  on the edge is:

.. math::

  \int_{\Gamma} \phi_i \left(n_f \cdot n_s - (cos\ \theta_s - c_T Ca_L)\right) d \Gamma = 0

  

|

  where φ\ :sub:`i` is the shape function associated with node *i*.



--------------
**References**
--------------

No References.