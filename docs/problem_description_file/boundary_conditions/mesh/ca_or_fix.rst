*************
**CA_OR_FIX**
*************

::

	BC = CA_OR_FIX NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This boundary condition card allows the application of Gibb’s inequality condition in
conjuction with a contact angle. This allows for a point to be specified at which a
contact line will attach itself and no longer move. Up to that point, the contact line will
advance or recede with a specified fixed contact angle. When the contact line attaches,
its contact angle is allowed to vary permitting the user to include discontinuities in
surface slope as features of the problem. The Gibb’s condition also permits the contact
line to detach from its fixed point if the contact angle enters a certain range after
attaching. This boundary condition is applicable only to two-dimensional problems;
see *CA_EDGE_OR_FIX* for details on three dimensional implementations.

The <float_list> has seven values, with definition of the input parameters as follows:

=================== ===================================================================
**CA_OR_FIX**       Name of the boundary condition.
**NS**              Type of boundary condition (<bc_type>), where **NS** denotes
                    node set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (node set in
                    EXODUS II) in the problem domain.
<float1>            θ\ :sub:`dcl`, dynamic contact angle, in radians.
<float2>            n\ :sub:`x`, x-component of outward-pointing wall surface normal.
<float3>            n\ :sub:`y`, y-component of outward-pointing wall surface normal.
<float4>            n\ :sub:`z`, z-component of outward-pointing wall surface normal.
<float5>            x\ :sub:`0`, x-coordinate of the point or feature at which the
					meniscus will pin.
<float6>            y\ :sub:`0`, y-coordinate of the point or feature at which the
					meniscus will pin.
<float7>            z\ :sub:`0`, z-coordinate of the point or feature at which the
					meniscus will pin.
=================== ===================================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = CA_OR_FIX NS 100 1.3 0. 1. 0. -0.5 1. 0.

-------------------------
**Technical Discussion**
-------------------------

The Gibb’s inequality condition is illustrated in the accompanying figure. The fixed
point is indicated by the plane, x = x\ :sub:`0`. Initially, the contact line is far from this point as
the condition at the contact line fixes the contact angle to the value θ\ :sub:`dcl`. However,
when the contact line approaches to within ε (1.e-6) of the fixed point, it attaches there
and stops moving. The contact angle condition is no longer enforced and the angle of
the free surface with respect to the solid normal vector is allowed to vary freely. The
other part of the Gibb’s inequality is illustrated (above) by the last sketch. Here, by
virture of the overall fluid mechanics, the contact angle withdraws until it is larger than
θ\ :sub:`dcl`. When this happens the contact line is no longer affixed at x = x\ :sub:`0` and is allowed to
move freely. Once again the contact angle condition is enforced.

.. figure:: /figures/056_goma_physics.png
	:align: center
	:width: 90%

	Contact angles and Gibb’s inequality condition in Goma, for the
	special case when the meniscus is moving along a surface of constant x.

Also, please see the important note under the BC = CA card regarding the convention
used for specifying wall and free surface normal vectors.

----------
**Theory**
----------

The principle behind this condition applies when a contact line encounters a sharp
feature on a surface. The feature from a distance might appear as a sharp corner at
which the meniscus/contact line prefers to locate rather than undergo continued
migration. Actually on a smaller scale, the corner feature is not infinitely small, and the
contact line undergoes no perceptable movement on the macroscale in order to satisfy a
true contact angle. Rather than resolving this feature with a fine mesh, it is an expedient
to pin the contact line there and allow it to take on any macroscale contact angle within
a certain range. The line can release again if the meniscus pulls the contact line
sufficiently to overcome specified bounds.


--------------
**References**
--------------

No References.