*****************
**NORM_FORCE_RS**
*****************

::

	BC = NORM_FORCE_RS SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR REALSOLID)**

This boundary condition card applies a force per unit area (traction) on a real-solid in a
*TOTAL_ALE* mesh region (see *Mesh Motion* card). The force per unit area is applied
uniformly over the boundary delineated by the side set ID. The applied traction is of
course a vector. Unlike the *FORCE_RS* boundary condition card, the vector traction
here is defined in normal-tangent vector basis. Definitions of the input parameters are
as follows:

================= =============================================================
**NORM_FORCE_RS** Name of the boundary condition (<bc_name>).
**SS**            Type of boundary condition (<bc_type>), where **SS**
                  denotes side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set
                  in EXODUS II) in the problem domain.
<float1>          Normal component of traction in units of force/area.
<float2>          Tangential component of traction in units of force/area.
<float3>          Second tangential component of traction in units of
                  force/area (in 3-D).
================= =============================================================

This card actually applies a traction that is then naturally integrated over the entire side
set of elements. Hence, the units on the floating point input must be force/area.

------------
**Examples**
------------

The following is a sample input card:
::

     BC = NORM_FORCE_RS SS 10 0. 1.0 1.0

This card results in a vector traction to the real-solid in a *TOTAL_ALE* mesh motion
type (not the mesh) defined by
:math:`\vec F` = 0.0(:math:`\vec n`) + 1.0(:math:`\vec t_1`) + 1.0(:math:`\vec t_2`) applied to the side set
boundary delineated by flag 10. The normal vector is defined as the outward pointing
normal to the surface. For internal surfaces defined by side sets which include both
sides of the interface, this condition will result in exactly a zero traction, i.e., internal
surface side sets must be attached to one element block only to get a net effect.

-------------------------
**Technical Discussion**
-------------------------

It is important to note that this boundary condition can only be applied to *TOTAL_ALE*
mesh motion types (cf. *Mesh Motion* card). As an example of how this boundary
condition card is used, consider the need to apply some load pressure uniformly on the
inside of a solid-membrane (like a pressurized balloon). In more advanced usage, one
could tie this force to an augmenting condition on the pressure, as dictated by the ideal
gas law.

This boundary condition is not used as often as the *FORCE_RS* or *FORCE_USER_RS*
counterparts.



--------------
**References**
--------------

No References.