***************
**SHARP_CA_2D**
***************

::

	BC = SHARP_CA_2D SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to impose a contact angle on a boundary when using
*Level Set Interface Tracking*. It can only be used for two-dimensional problems.

A description of the input parameters follows:

=========== ===============================================================
**FILL_CA** Name of the boundary condition.
**SS**      Type of boundary condition (<bc_type>), where **SS** denotes
            side set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (side set in
            EXODUS II) in the problem domain.
<float>     :math:`\theta`, the contact angle imposed, in degrees.
=========== ===============================================================

------------
**Examples**
------------

An example:
::

   BC = SHARP_CA_2D SS 10 30.0

-------------------------
**Technical Discussion**
-------------------------

This boundary condition must be used in conjunction with the *VELO_SLIP_FILL* or
*VELO_SLIP_LS* boundary condition. These latter conditions permits the fluid to slip in
the vicinity of the contact line. The SHARP_CA_2D acts by imposing a force on the
momentum equation. The size of this force is more or less in proportion between the
actual contact angle on the boundary and the value specified on the card and scales
directly with the applied surface tension material parameter. In this manner, it is very
similar to the FILL_CA boundary condition.

The manner in which is applied differs. In this case, the applied force is not distributed
around the contact line using a smooth delta function weighting in a weak integrated
context, but instead the delta function is used to resolve the line integral and the 
force is applied directly at a point on the sideset set. Hence, this boundary condition 
is most appropriate for use in conjunction with subelement integration which performs a
similar transformation of the volumetric surface tension source terms. Further, the
logic use to identify the point of application on the boundary functions only in 
twodimensions. Hence, this boundary condition is stricly limited to two-dimensional
problems.



--------------
**References**
--------------

No References. 