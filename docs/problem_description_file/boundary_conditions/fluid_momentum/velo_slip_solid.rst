*******************
**VELO_SLIP_SOLID**
*******************

::

	BC = VELO_SLIP_SOLID SS <bc_id> <integer_list> <float1> [integer3, float2]

-----------------------
**Description / Usage**
-----------------------

**(WIC/ROTATED MOMENTUM)**

This boundary condition is similar in function to the *VELO_SLIP* condition in that it
permits a tangential velocity in a fluid phase to be proportional to the shear stress at the
boundary. This boundary condition allows for this type of slip to occur at the interface
between a fluid material and a *LAGRANGIAN* or *TALE* solid material. The velocity of
the solid substrate is obtained automatically from the motion of the solid material,
including advection of the stress-free state. As in the case of the *VELO_SLIP* condition,
this condition also permits the user to vary the slip coefficient depending upon the
distance from a specified point in the mesh. The variable slip model can only be used in
two-dimensional problems.

The <integer_list> has two values; the definitions of the input parameters and their
significance in the boundary condition parameterization is described below:

=================== =================================================================
**VELO_SLIP_SOLID** Name of the boundary condition.
**SS**              Type of boundary condition (<bc_type>), where **SS**
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set
                    in EXODUS II) in the problem domain. This should be
                    an internal sideset defined at the interface between solid
                    and liquid material blocks.
<integer1>          The element block id defining the solid material phase.
<integer2>          The element block id defining the liquid material phase.
<float1>            :math:`\beta`, the slip coefficient. The inverse of 
                    :math:`\beta` defines the
                    scaling between stress and slip. Hence, for small values
                    of :math:`\beta`, large shear stresses are needed for a 
                    given amount of slip, and conversely, for large values of 
                    :math:`\beta`, the amount
                    of stress needed for the same degree of slip decreases
                    (see below for a more rigorous description).
[integer3]          :math:`N_{cl}`, a single-node node set identification number.
                    When the variable coefficient slip relation is used, distance 
                    is measured relative to this node (see discussion below).
                    Normally, this node set represents the location of the
                    dynamic contact line. Note that this option is generally
                    only used in two-dimensional simulations.
[float2]            :math:`\alpha`, the distance scale in the variable slip model 
                    (see the discussion below). Both :math:`N_{cl}` and 
                    :math:`\alpha` should be present to
                    activate the variable slip model.
=================== =================================================================

------------
**Examples**
------------

The following is a sample card:
::

     BC = VELO_SLIP_SOLID SS 20    2 1 0.001 0.0 4 0.01

This boundary condition sets the slip coefficient between solid material 2 and liquid
material 1 to be 0.001 except in the vicinity of the nodeset 4 (a single node) where the
variable model is used.

-------------------------
**Technical Discussion**
-------------------------

* The general form of this boundary condition is

.. figure:: /figures/093_goma_physics.png
	:align: center
	:width: 90%

where :math:`\tau` is the deviatoric portion of the fluid stress tensor,  :math:`\beta`
is the Navier slip
coefficient and :math:`v_{sfs}` is the velocity of the solid surface stress-free state, with :math:`F_m` the
deformation gradient tensor; this motion includes any rigid solid body motion and
any superimposed deformation velocity.

* It is worthwhile noting that, unlike the *VELO_SLIP* condition, this condition is
  actually a rotated condition. It is applied to the tangential component of the rotated
  momentum equations weakly. This means that the normal component of the
  momentum equation is not affected by this boundary condition. Normally, some
  sort of no-penetration condition must accompany this boundary condition for this
  reason.

* The reader is referred to the documentation of the variable slip coefficient model to
  apply slip near contact lines under the *VELO_SLIP* boundary condition.




.. TODO - Line 78 contains a photo that needs to be exchanged for the equation.