***********
**LS_CA_H**
***********

::

	BC = LS_CA_H SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/SCALAR CURVATURE)**

This boundary condition is used only in conjunction with level set interface tracking
and the LS_CAP_CURVE embedded surface tension source term. Its function is
impose a contact angle condition on that boundary.

A description of the input parameters follows:

=========== ====================================================================
**LS_CA_H** Name of the boundary condition.
**SS**      Type of boundary condition (<bc_type>), where **SS** denotes
            side set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (side set in
            EXODUS II) in the problem domain.
<float>     A float value that is the imposed contact angle in degrees.
=========== ====================================================================

------------
**Examples**
------------

An example:
::

   BC = LS_CA_H SS 10 45.0

-------------------------
**Technical Discussion**
-------------------------

The projection equation operator for solving for the curvature degree of freedom from a
level set field is a Laplacian. It is standard to integrate these operators by parts but in
the process one always generates a boundary integral. In this case the integral takes the
form:

.. figure:: /figures/197_goma_physics.png
	:align: center
	:width: 90%

where :math:`n_w` is the wall surface normal and :math:`n_{fs}` is the normal to free surface (zero contour
of the level set function ). This is a convenient event because it allows us to impose a
contact angle condition on a sideset using this boundary integral by making the
assignment

.. figure:: /figures/198_goma_physics.png
	:align: center
	:width: 90%

where :math:`\theta` is the contact angle specified on the card.

The effect of this boundary condition is impose a disturbance in the curvature field
near the boundary that has the effect of accelerating or decelerating the fluid near the
wall in response to whether the actual contact angle is greater or less than the imposed
value. Thus, over time, given no other outside influences, the contact angle should
evolve from its initial value (that presumably is different than the imposed value) to the
value imposed on this card. The user should expect that the contact angle will
instantaneously jumped to the imposed value.



--------------
**References**
--------------

No References. 

..
	TODO - Lines 49 and 58 have pictures that need to be exhanged with equations.



