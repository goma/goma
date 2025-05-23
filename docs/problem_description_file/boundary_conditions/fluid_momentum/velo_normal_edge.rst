********************
**VELO_NORMAL_EDGE**
********************

::

	BC = VELO_NORMAL_EDGE SS <bc_id1> <bc_id2> <float>

-----------------------
**Description / Usage**
-----------------------

**(PCC-EDGE/ROTATED MOMENTUM)**

This boundary condition card is used to specify the normal velocity component on a
dynamic contact line in three-dimensions. The velocity component is normal to the
contact line in the plane of the web and is equal to :math:`V_n`. The free-surface side set should
always be <bc_id1>, the primary side set, and the web side set should be <bc_id2>, the
secondary side set. Usually, this boundary condition is used to model dynamic contact
lines in three dimensions and is usually found in conjunction with a
*VELO_TANGENT_EDGE* card, a *VAR_CA_EDGE* or *CA_EDGE* card as explained
below.

Definitions of the input parameters are as follows:

===================== ===========================================================
**VELO_NORMAL_EDGE**  Name of the boundary condition.
**SS**                Type of boundary condition (<bc_type>), where **SS**
                      denotes side set in the EXODUS II database.
<bc_id1>              The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain. This is the
                      primary side set defining the edge and should also be
                      associated with the capillary free surface if used in the
                      context of a dynamic contact line.
<bc_id2>              The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (side set
                      in EXODUS II) in the problem domain. Together with
                      <bc_id1>, this secondary side set defines the edge/curve
                      on which the boundary condition applies as the
                      intersection of the two side sets. In problems involving
                      dynamic contact lines, this side set should correspond to
                      the moving substrate.
<float>               :math:`V_n`, a parameter supplying the imposed normal velocity
                      component. This component is taken normal to the edge
                      curve parallel to <bc_id2>. See below for a more
                      detailed description.
===================== ===========================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_NORMAL_EDGE SS 5 4 0.0

This card sets the normal-to-contact line component of the velocity to zero along the
curve defined by the intersections of side set 5 and 4.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition imposes a point
  collocated constraint of the form:

.. math::

  n_{cl} \cdot \left(v - v_m\right) = V_n

  

where *v* is the fluid velocity, :math:`v_m` is the mesh
velocity and :math:`n_cl` is the normal to the contact
line in the plane of <bc_id2>. The sketch at
right depicts the orientation of this latter
vector. Note that the collocation points for this boundary condition only are not the
nodes on the edge curve but integration points in each of the edge elements. The
reason for this is historical and uninteresting from a user point of view.

.. figure:: /figures/080_goma_physics.png
	:align: center
	:width: 90%

* This boundary condition is used almost exclusive in problems involving dynamic
  contact lines in three dimensions. Imposition of wetting line physics is a 
  difficult
  problem in modeling situations involving dynamic contact lines. In twodimensions,
  the assumption is often made that the effect of any wetting line force
  is to locally produce a condition in which the fluid velocity at the contact line is
  zero *in the laboratory reference frame*. That is to say, that at the contact line noslip
  between fluid and moving substrate is not enforced and instead a zero velocity
  condition is imposed. In this way, the difficult-to-model wetting line forces are not
  included directly, but instead are included by their effect on the velocity. One
  might argue with this model, and many do, but as a practical approach, this has
  been shown to work well.

  Generalizing this notion into three dimensions is the primary motivation for this
  boundary condition. In the case of a dynamic contact line that is a curve in three
  dimensions, it is not correct to simply set all velocity components to zero because
  that would imply that the wetting forces act equally in all three directions. It is
  more reasonable to say that the wetting forces can act only in a direction normal to
  the contact line in the plane of the substrate. Therefore, the correct generalization
  of the wetting line model described in the previous paragraph is to set the velocity
  component normal to the contact line in the plane of the substrate to zero. This is
  done by using the *VELO_NORMAL_EDGE* boundary condition with :math:`V_n` set to
  zero. In the case of a transient problem, it is necessary to add the qualifier,
  “relative to the mesh motion.” This accounts for the mesh motion velocity in the
  constraint equation. See Baer, et.al. (2000) for a more complete discussion of this
  wetting line model.

* Generally, a *VELO_NORMAL_EDGE* card must be accompanied by other
  boundary conditions for a correct application. Firstly, since
  *VELO_NORMAL_EDGE* forces the velocity vector to be parallel to the contact
  line (at least in steady state), the *KINEMATIC* condition on any free surface
  attached to the contact line will overspecify the problem at the contact line. For
  this reason, it is generally the case that a *CA_EDGE, VAR_CA_EDGE* or
  *VAR_CA_USER* (or their variants) should also be present for the contact line.
  These boundary conditions replace the *KINEMATIC* card on the mesh at the
  contact line.

  In addition, a *VELO_TANGENT_EDGE* card should be present to enforce no-slip
  between fluid and substrate *in the tangential direction*. Also it should be
  recognized that *VELO_NORMAL_EDGE* will not override other Dirichlet
  conditions on the substrate side set. Typically, the latter are used to apply no slip
  between fluid and substrate. If such conditions are used over the entirety of the
  substrate side set, both *VELO_NORMAL_EDGE* and *VELO_TANGENT_EDGE*
  conditions applied at the contact will be discarded.

  There are two potential solutions to this. First, the substrate region could be
  divided into two side sets, a narrow band of elements adjacent to the contact line
  and the remainder of substrate region. In the narrow band of elements, the no slip
  condition is replaced by a *VELO_SLIP* card with the substrate velocity as
  parameters. This allows the velocity field to relax over a finite region from the
  velocity imposed at the contact line to the substrate field. The second method uses
  only a single side set for the substrate region, but replaces the Dirichlet no slip
  boundary conditions with a penalized *VELO_SLIP* condition. That is, the slip
  parameter is set to a small value so that no slip is effectively enforced, but within
  the context of a weakly integrated condition. Since the *VELO_NORMAL_EDGE*
  and *VELO_TANGENT_EDGE* cards are strongly enforced on the contact lines, the
  *VELO_SLIP* card will be overridden in those locations and the velocity field will
  deviate appropriately from the substrate velocity.



--------------
**References**
--------------

Baer, T.A., R.A. Cairncross, P.R.Schunk, R.R. Rao, and P.A. Sackinger, “A finite
element method for free surface flows of incompressible fluids in three dimensions.
Part II. Dynamic wetting lines.” IJNMF, 33, 405-427, (2000).