**************
**DXYZDISTNG**
**************

::

	BC = {DXDISTNG | DYDISTNG | DZDISTNG} SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(PCC/MESH)**

This boundary condition card is used to specify a distinguishing condition for mesh
motion based on an isotherm, viz. the distinguishing condition forces the mesh
boundary to which it is applied to take on a position such that the temperature is
constant and at the specified value, all along the boundary. Although of the same
mathematical form as the *DISTNG* boundary condition, this condition does not force a
boundary rotation of the vector mesh residuals. Instead, it is recommended that the
condition be chosen such that the predominant direction of the normal vector is close to
one of the three Cartesian coordinates, X, Y, or Z. For example, if the boundary in
question is basically oriented so that the normal vector is mostly in the positive or
negative Y-direction, then *DYDISTNG* should be chosen. Definitions of the input
parameters are as follows:

==================================== =============================================================
**{DXDISTNG | DYDISTNG | DZDISTNG}** Eight-character boundary condition name (<bc_name>)
                                     that defines the distinguishing condition, where:
                                     
                                     	* **DXDISTNG** - X condition
                                     	* **DYDISTNG** - Y condition
                                     	* **DZDISTNG** - Z condition
**SS**                                Type of boundary condition (<bc_type>), where **SS**
                                      denotes side set in the EXODUS II database.
<bc_id>                               The boundary flag identifier, an integer associated with
                                      <bc_type> that identifies the boundary location (side set
                                      in EXODUS II) in the problem domain.
<float>                               Value of temperature isotherm. If one wanted to apply a
                                      variable temperature, e.g. as a function of the
                                      concentration, it is suggested that the user-defined
                                      boundary conditions be used.
==================================== =============================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = DYDISTNG SS 123 273.0

This card forces the boundary defined by EXODUS II side set number 123 to conform
to the isotherm temperature of 273.0. Most importantly, the y-component of the mesh
equation residuals is replaced by this condition.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form of this distinguishing condition is as follows:

.. math::

   T - T_{\mathrm{mp}} = 0

where :math:`T_{\mathrm{mp}}` is the specified temperature parameter. This condition has been used
extensively for macroscale and microscale melting problems, whereby one needs to
distinguish a molten region from a solidified or mushy region with liquidus and solidus
temperatures. In three dimensions usage needs to be completed with a companion *ROT*
input card which directs the equation application of the condition, even though
rotations are not actually performed.

As a bit of software trivia, this is the first distinguishing condition ever written in
*Goma*, and one of the first boundary conditions, period.




.. 
	TODO - The image on line 64 needs to be replaced by the actual equation.