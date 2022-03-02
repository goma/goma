***********
**FILL_CA**
***********

::

	BC = FILL_CA SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition is used to impose a contact angle on a boundary when using
*Level Set Interface Tracking*.

A description of the input parameters follows:

============ ===============================================================
**FILL_CA**  Name of the boundary condition.
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<float>      :math:`\theta`, the contact angle imposed, in degrees.
============ ===============================================================

------------
**Examples**
------------

An example:
::

    BC = FILL_CA SS 10 30.0

-------------------------
**Technical Discussion**
-------------------------

This boundary condition must be used in conjunction with the *VELO_SLIP_FILL*
boundary condition. This latter condition permits the fluid to slip in the vicinity of the
contact line. The *FILL_CA* acts by imposing a force on the momentum equation. The
size of this force is more or less in proportion between the actual contact angle on the
boundary and the value specified on the card. This force is applied as a weakly
integrated condition and if the *VELO_SLIP_FILL* condition is not present, the
*FILL_CA* will be overwritten and *ipso facto* absent.



--------------
**References**
--------------

No References.