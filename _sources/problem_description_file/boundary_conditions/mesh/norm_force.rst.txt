**************
**NORM_FORCE**
**************

::

	BC = NORM_FORCE SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MESH)**

This boundary condition card applies a force per unit area (traction) on a Lagrangian
mesh region. The force per unit area is applied uniformly over the boundary delineated
by the side set ID. The applied traction is of course a vector. Unlike the *FORCE*
boundary condition card, the vector traction here is defined in normal-tangent vector
basis. Definitions of the input parameters are as follows:

=============== ==================================================================
**NORM_FORCE**  Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set.
<bc_id>         The boundary flag identifier, or a side set number which
                is an integer that identifies the boundary location (side
                set in EXODUS II) in the problem domain.
<float1>        Normal component of traction in units of force/area.
<float2>        Tangential component of traction in units of force/area.
<float3>        Second tangential component of traction in units of
                force/area (in 3-D).
=============== ==================================================================

This card actually applies a traction that is then naturally integrated over the entire 
side set of elements. Hence, the units on the floating point input must be force/area.

------------
**Examples**
------------

Following is a sample card:
::

     BC = NORM_FORCE SS 10   0.   1.0   1.0

This card results in a vector traction defined by **EQUATION** being
applied to the side set boundary delineated by the number 10. The normal vector is
defined as the outward pointing normal to the surface. For internal surfaces defined by
side sets which include both sides of the interface, this condition will result in exactly 
a zero traction, i.e., internal surface side sets must be attached to one element block 
only to get a net effect.

-------------------------
**Technical Discussion**
-------------------------

**Important note**: this boundary condition can only be applied to *LAGRANGIAN,
DYNAMIC_LAGRANGIAN* or *ARBITRARY* mesh motion types (cf. *Mesh Motion* card).
For real-solid mesh motion types, refer to *NORM_FORCE_RS*. Furthermore, it is rare
and unlikely that this boundary condition be applied to *ARBITRARY* mesh motion
regions. An example application of this boundary condition card is to apply some load
pressure uniformly on the inside of a solid-membrane (like a pressurized balloon). In
more advanced usage, one could tie this force to an augmenting condition on the
pressure, as dictated by the ideal gas law.

This boundary condition is not used as often as the *FORCE* or *FORCE_USER*
counterparts.




..
	 TODO - Where it says "**EQUATION**" in line 46 there is supposed to be an equation that needs to be written.