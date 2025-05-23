****************
**VELO_TANGENT**
****************

::

	BC = VELO_TANGENT SS <bc_id> <integer> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MOMENTUM)**

This boundary condition is used to specify strongly the component of velocity
tangential to the side set. An added feature is the ability to relax the condition near a
point node set according to supplied length scale and slipping parameters. This has
application to problems involving moving contact lines. Note that this boundary
condition is applicable only to two-dimensional problems and will result in an error if it
is used in a three-dimensional context.

The <float_list> has three parameters; definitions for all input parameters is as follows:

================ =====================================================================
**VELO_TANGENT** Name of the boundary condition.
**SS**           Type of boundary condition (<bc_type>), where **SS**
                 denotes side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set
                 in EXODUS II) in the problem domain.
<integer>        :math:`N_{cl}`, parameter that identifies a single-node node
                 set that coincides with the location in the model of the 
                 moving contact line. Distances in the slipping model are
                 computed relative to the location of this node. When the
                 slipping model is not used, this parameter can safely be
                 set to zero. Another toggle setting can be triggered by
                 setting this integer to -1; with this the *VELO_TANGENT*
                 condition is kept at a rolling motion dynamic contact
                 line. (See FAQ below on rolling motion conditions.)
<float1>         :math:`v_t`, a parameter specifying the value of the tangent
                 velocity component. The component direction is *n* × *k*
                 where *k* is the z-component unit vector.
<float2>         :math:`\beta`, a parameter specifying the coefficient for
                 slip velocity (see model below); setting :math:`\beta`
                 to zero disables the slipping model.
<float3>         :math:`\alpha`, a parameter specifying the length scale for the
                 position dependent slip (see model below); setting 
                 :math:`\alpha` to zero disables the slipping model.
================ =====================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VELO_TANGENT SS 10 100 0.0 1.0 0.1

-------------------------
**Technical Discussion**
-------------------------

* Most often this boundary condition is used only to set the tangential speed on a
  side set because simpler Dirichlet conditions are not appropriate. An example is a
  sloping fully-developed inlet plane which does coincide with a coordinate axis. In
  this case, this boundary condition would be used to set the tangential velocity to be
  zero. The constraint applied at node *i* is as follows:

.. math::

  \int_\Gamma \phi_i \left(t \cdot v - v_t \right)d \Gamma = 0



* Alternatively, a dynamic contact line might be present in the problem and it is
  desirable that this condition be relaxed near the position of this contact line. This
  can be done by supplying non-zero values for :math:`\alpha` and :math:`\beta`. In this case, the constraint
  that is applied at the :math:`i^{th}` node on the boundary is:

.. math::

  \int_\Gamma \phi_i \left(t \cdot v - \beta \dot x e^{-\alpha d} - v_t\right)d \Gamma = 0

  

|

  in which *d* is the straightline distance to the node attached to <:math:`N_{cl}`> and :math:`\dot{x}` is the
  velocity vector of the mesh. It should be recognized that for steady state problems
  the mesh motion is by definition always zero so this constraint reverts to the
  previous expression.


--------
**FAQs**
--------

**Rolling Motion Conditions for high Capillary number dynamic wetting**. Often times it
is desirable to model a case of dynamic wetting for which the conditions result in a high
capillary number. At this limit, it is well known that a contact angle specification is in
fact an overspecification. *Goma* has always been able to model this case, except
recently some changes have been made to allow for the combination of conditions at a
dynamic contact line to be controlled. It should be stressed that all finite capillary
number cases still work as always. This FAQ addresses the special case in which you
desire to specify no-slip right up to the contact line. In most cases a *VELO_SLIP* card
or outright setting the velocity components to zero at the moving contact line in order
to impart slip will circumvent the issue taken up here.

The figure below diagrams the situation:

.. figure:: /figures/085_goma_physics.png
	:align: center
	:width: 90%

Basically the web in this example corresponds to side set 5 and the free surface to side
set 4. The conditions we desire in the vicinity of the contact line are as follows:
::

	$web surface
    BC = VELO_TANGENT SS 5   0   {web_sp}   0.0   0.0
    BC = VELO_NORMAL SS 5   0.0
    BC = GD_PARAB   SS 5   R_MESH2   0   MESH_POSITION1   0   0. 0. 1.
    BC = GD_PARAB   SS 5   R_MESH2   0 MESH_POSITION2   0 0.   {2*roll_rad} 1.
    $   upstream heel
    BC = KINEMATIC SS 4 0.
    BC = CAPILLARY   SS   4 {inv_cap}   0.0   0.0

Notice how there is no contact angle specified and even with the CAPILLARY card, the
effect of,

VELO_NORMAL, surface tension is very small. The desired set of conditions that
should be applied at the dynamic contact line are as follows:
::

	At node 1:
    R_MOMENTUM1   gets VELO_NORMAL   from SS 5, CAPILLARY   from SS 4,
    R_MOMENTUM2   gets VELO_TANGENT   from SS 5, CAPILLARY   from SS 4,
    R_MESH1   gets KINEMATIC   from SS 4,
    R_MESH2   gets GD_PARAB   from SS 5, GD_PARAB   from SS 5,

This clearly shows that at the contact line, which happens to be node number 1 as
shown by this clip from the BCdup.txt file resulting from the run, both
*VELO_NORMAL* and *VELO_TANGENT* cards are applied, which implies no-slip. This
is the so-called rolling-motion case (or tank-tread on a moving surface) in which the
“kinematic paradox” is no longer a paradox. That is, both the *KINEMATIC* condition
on the free surface and the no-slip condition on the substrate can be satisfied without
loss or gain of mass through the free surface (see Kistler and Scriven, 1983). In order to
make sure that both the combination above is applied, a “-1” must be placed in the first
integer input of the *VELO_TANGENT* card, vis.,
::

	BC = VELO_TANGENT SS 5   -1 {web_sp}   0.0   0.0

This integer input slot is actually reserved for a variable slip coefficient model and is
normally used to designate the nodal bc ID of the contact line. In this case of no-slip, it
is not needed so we added this special control. If the following card is issued:
::

	BC = VELO_TANGENT SS 5   0 {web_sp}   0.0   0.0

then the following combination results:
::

	At node 1:
    R_MOMENTUM1 gets   VELO_NORMAL   from SS 5, CAPILLARY   from SS 4,
    R_MOMENTUM2 gets   CAPILLARY   from SS 4,
    R_MESH1 gets   KINEMATIC   from SS 4,
    R_MESH2 gets   GD_PARAB   from SS 5, GD_PARAB   from SS 5,

which is desired in the case for which a contact angle and liquid slip is applied.

--------------
**References**
--------------

Kistler, S. F. and Scriven, L. E. 1983. Coating Flows. In Computational Analysis of
Polymer Processing. Eds. J. A. Pearson and S. M. Richardson, Applied Science
Publishers, London.

.. TODO - In lines 70 and 79 the equations need to replace the pictures.