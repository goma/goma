**********
**DISTNG**
**********

::

	BC = DISTNG SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This boundary condition card is used to specify a distinguishing condition for mesh
motion based on an isotherm, viz. the distinguishing condition forces the mesh
boundary to which it is applied to take on a position such that the temperature is
constant and at the specified value, all along the boundary. This condition causes the 
vector mesh motion equations (viz. *mesh1, mesh2, and mesh3* on *EQ* cards) to be
rotated into normal-tangential form. In two dimensions, this condition is applied to the
normal component automatically; in three dimensions it is suggested to put it on the
normal component, as specified by the *ROT* conditions. Definitions of the input
parameters are as follows:

================ =============================================================================
**DISTNG**       Name of the boundary condition (<bcname>).
**SS**           Type of boundary condition (<bc_type>), where SS denotes
                 side set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (side set in 
                 EXODUS II) in the problem domain.
<float>          Value of temperature/isotherm. To apply a variable
                 temperature, e.g., as a function of the concentration, it is
                 suggested that the user-defined boundary conditions be
                 used, like *SPLINE* or *GEOM*.
================ =============================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = DISTNG SS 123 273.0

This card forces the boundary defined by EXODUS II side set number 123 to conform
to the isotherm temperature of 273.0.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form of this distinguishing condition is as follows:

.. math::

   T - T_{\mathrm{mp}} = 0

where :math:`T_{\mathrm{mp}}` is the specified temperature parameter. This condition has been used
extensively for macroscale and microscale melting problems, whereby one needs to
distinguish a molten region from a solidified or mushy region with liquidus and solidus
temperatures. In three dimensions, usage needs to be completed with a companion *ROT*
input card which directs the equation application of the condition.


--------
**FAQs**
--------

Continuation Strategies for Free Surface Flows In free surface problems, there exists
one or more boundaries or internal surfaces whose position(s) are unknown a *priori*. As
such, the geometry of the problem becomes part of the problem and must be
determined together with the internal physics. Most problems of this sort cannot be
solved with a trivial initial guess to the solution vector, mainly because the conditions
which determine the surface position are closely coupled to the active physics in the
bulk. Thus, these problems require continuation (zero or higher order) to achieve a
converged solution to a desired state. The continuation strategy typically involves
turning on and off the conditions which distinguish the position of the free surface(s);
one such strategy is described in this FAQ.

Distinguishing conditions in *Goma* serve two purposes: (1) they can be used to locate a
surface whose position depends on internal and interfacial transport phenomena, and
(2) they can be used to prescribe solid boundary position or motion. The first type of
condition contains field variables needed to locate the interface or free surface position,
and hence ties the mesh motion to the problem physics, i.e., mass, momentum, and
energy transport phenomena. Currently, the side-set boundary conditions of type
*DISTNG, KINEMATIC,* and *KIN_LEAK* fall into this class. The second type of
condition requires only geometrical information from the mesh, and, although
geometrically couples the mesh motion to the problem physics, it tends not to be so
tightly coupled. Currently, boundary conditions *PLANE, PLANEX, PLANEY, PLANEZ,
SPLINE, SPLINEX, SPLINEY,* and *SPLINEZ* fall into this class.

In two dimensions, there is no need to use *PLANEX, PLANEY, PLANEZ, SPLINEX,
SPLINEY,* and *SPLINEZ*. Because the code automatically rotates the mesh residual
equations and the corresponding Jacobian entries into normal-tangential form on the
boundary, *SPLINE, PLANE,* and *DISTNG* are the only cards required to specify the
position of the boundary. Currently, in three dimensions, the logic for the same rotation
concept is not totally functional, and one must use the *PLANEX*, etc. cards to designate
which component of the mesh stress residual equation receives the distinguishing
conditions.

If cards *DISTNG, KINEMATIC,* and *KIN_LEAK*, i.e., distinguishing conditions of type
1, are absent in any simulation, then any initial guess for the transport field equations,
i.e., energy and momentum, has a chance of converging, as long as the initial mesh
displacement guess is within the radius of convergence of the mesh equations and
associated boundary conditions. For example, if the side sets of the EXODUS II
database mesh correspond somewhat closely to what is prescribed with *PLANE* and
*SPLINE*-type conditions, then an initial guess of the NULL vector has a good chance of
converging, so long as the velocities and temperatures are within “converging
distance.”

When conditions from the first class are present, i.e., either *DISTNG, KIN_LEAK* or
*KINEMATIC*, then the following procedure should be followed:

* Set the keyword for the *Initial Guess* character_string to **zero**, **one**, or **random**.

* Obtain a solution (run *Goma*) with the initial guess for the free surfaces
  distinguished as *KINEMATIC* (or other) coming from the EXODUS II database,
  but without the *KINEMATIC* (or other) card(s). That is, “fix” those surfaces with
  either a *PLANE* or *SPLINE* command, or simply place no distinguishing condition
  on them (this works only if the grid has not been previously “stressed”, i.e., all the
  displacements are zero). The rest of the “desired” physics should be maintained. If
  any surface is distinguished as *KINEMATIC*, then it is highly advantageous to
  place a *VELO_NORMAL* condition on that surface for startup, and set the
  corresponding floating point datum to zero. This effectively allows the fluid to
  “slip” along that boundary as if it were a shear free condition.

* Set the keyword in the *Initial Guess* character_string to **read**.

* Copy the file named in *SOLN file* into the file named in *GUESS file*.

* Release the free boundaries by taking off any current distinguishing condition
  cards and adding the appropriate *KINEMATIC* (or other) card. Adjust all other
  boundary conditions appropriately.

* Run *Goma*, using a relaxed Newton approach (factor less than unity but greater
  than zero - e.g., 0.1) for complex flows.

When dealing with material surface boundaries distinguished by the kinematic
boundary condition, the nature of that condition requires a non-zero and substantial
component of velocity tangent to the surface upon start-up. In this case, it can be
advantageous to use the *VELO_TANGENT* card to set the velocity along the free
surface to some appropriate value prior to releasing the free surface (in the third step
above). Of course this card will be removed in subsequent steps. Also, although not
necessary, a smooth, “kinkless”, initial guess to the free surface shape is helpful
because it reduces the amount of relaxation required on the Newton iteration.

Obtaining start-up solutions of most coating flow configurations is still an art. The best
way to start up a coating flow analysis may be to acquire a “template” developed from
a previous analysis of some closely related flows.

--------------
**References**
--------------

Allen Roach’s or Randy’s ESR tutorials. Perhaps these need to be put into the
repository.

.. 
	TODO - The picture in line 56 needs to be replaces by the actual equation.