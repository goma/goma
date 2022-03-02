************
**FORCE_RS**
************

::

	BC = FORCE_RS SS <bc_id> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR REALSOLID)**

This boundary condition card applies a force per unit area (traction) on a real-solid
material region (as opposed to a Lagrangian solid region), as is the case with
*TOTAL_ALE* mesh motion type (see *Mesh Motion* card). The force per unit area is
applied uniformly over the boundary delineated by the side set ID. The applied force is
of course a vector. Definitions of the input parameters are as follows:

============ ==============================================================
**FORCE_RS** Name of the boundary condition (<bc_name>)
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<float1>     X-component of traction in units of force/area.
<float2>     Y-component of traction in units of force/area.
<float3>     Z-component of traction in units of force/area.
============ ==============================================================

------------
**Examples**
------------

Following is a sample card:
::

     BC = FORCE_RS SS 10 0. 1.0 1.0

This card results in a vector traction defined by :math:`\vec F` = 0.0(:math:`\vec ex`) + 1.0(:math:`\vec ey`) + 1.0(:math:`\vec ez`) applied
to the side set boundary delineated by flag 10, where the element block bounded by this
boundary is of a *TOTAL_ALE* mesh motion type.

-------------------------
**Technical Discussion**
-------------------------

It is important to note that this boundary condition can only be applied to *TOTAL_ALE*
mesh motion types (cf. *Mesh Motion* card). (see *FORCE* for all other mesh motion
types). Furthermore, it is rare and unlikely that this boundary condition be applied to
*ARBITRARY* mesh motion regions. As an example of how this boundary condition card
is used, consider the need to apply some load pressure to a real solid of a *TOTAL_ALE*
region, like a rubber roller, so as to squeeze and drive flow in a liquid region. Some of
the usage tutorials cited below will direct you to some specifics.


--------
**FAQs**
--------

On internal two-sided side sets, this boundary condition results in double the force in
the same direction.

--------------
**References**
--------------

GT-005.3: THE NEW TOTAL-ARBITRARY-LAGRANGIAN-EULERIAN (TALE)
CAPABILITY and its applicability to coating with/on deformable media, August 6,
1999, P. R. Schunk

SAND2000-0807: TALE: An Arbitrary Lagrangian-Eulerian Approach to Fluid-
Structure Interaction Problems, P. R. Schunk (May 2000) 