*********
**FORCE**
*********

::

	BC = FORCE SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MESH)**

This boundary condition card applies a force per unit area (traction) on a Lagrangian
mesh region. The force per unit area is applied uniformly over the boundary delineated
by the side set ID. The applied force is of course a vector. Definitions of the input
parameters are as follows:

============ =======================================================================
**FORCE**    Name of the boundary condition (<bc_name>)
**SS**       Type of boundary condition (<bc_type>), where SS denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<float1>     X-component of traction in units of force/area.
<float2>     Y-component of traction in units of force/area.
<float3>     Z-component of traction in units of force/area.
============ =======================================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = FORCE SS 10 0. 1.0 1.0

This card results in a vector traction defined by **EQUATION** being
applied to the side set boundary delineated by the number 10.

-------------------------
**Technical Discussion**
-------------------------

**Important note**: this boundary condition can only be applied to *LAGRANGIAN,
DYNAMIC_LAGRANGIAN* or *ARBITRARY* mesh motion types (cf. *Mesh Motion* card).
For real-solid mesh motion types, refer to *FORCE_RS*. Furthermore, it is rare and
unlikely that this boundary condition be applied to *ARBITRARY* mesh motion regions.
An example application of this boundary condition card is to address the need to apply
some load pressure to a solid Lagrangian region, like a rubber roller, so as to squeeze
and drive flow in a liquid region.


--------
**FAQs**
--------

On internal two-sided side sets, this boundary condition results in double the force in
the same direction.

--------------
**References**
--------------

A MEMS Ejector for Printing Applications, A. Gooray, G. Roller, P. Galambos, K.
Zavadil, R. Givler, F. Peter and J. Crowley, Proceedings of the Society of Imaging
Science & Technology, Ft. Lauderdale FL, September 2001.

..
	 TODO - Where it says "**EQUATION**" in line 41 there is supposed to be an equation that needs to be written.