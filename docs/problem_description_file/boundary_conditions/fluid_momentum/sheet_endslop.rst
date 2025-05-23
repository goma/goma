******************
**SHEET_ENDSLOPE**
******************

::

	BC = SHEET_ENDSLOPE NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SPECIAL/VECTOR MOMENTUM)**

This boundary condition card is used to enforce a slope of a membrane surface (cf. to
be used in conjuction with *BC = TENSION_SHEET*) at its enpoints. . There are two
values to be input for the <float_list>; definitions of the input parameters are as
follows:

================== ========================================================
**SHEET_ENDFORCE** Name of the boundary condition (<bc_name>).
**NS**             Type of boundary condition (<bc_type>), where **NS**
                   denotes node set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (node
                   set in EXODUS II) in the problem domain.
<float1>           X-component of upstream idler point (see discussion
                   below)
<float2>           Y-component of upstream idler point.
================== ========================================================

------------
**Examples**
------------

The following is a sample input card using several APREPRO variables:
::

   BC = SHEET_ENDSLOPE NS 100   {sind(th2)}   {-cosd(th2)

This condition would enforce a slope equivalent to that defined between the
coordinates of the node at NS 100 with the point [sind(th2), -cosd(th2)].

-------------------------
**Technical Discussion**
-------------------------

Only two dimensional applications, viz. the nodeset is a single-node nodeset.

This is a single point nodeset boundary condition. Its function is to set the slope of the
web at the single point nodeset N. It does this by enforcing continuity of the slope of
the TENSION_SHEET sideset with the straight line that connects the point (X,Y) with
nodeset N. Thus, this boundary condition can be used to model the influence of an
upstream idler roller located at the point (X,Y). Indeed, this boundary condition has an
alternate name: IDLER_LOC.

This boundary condition exploits the natural boundary conditions associated with the
TENSION_SHEET formulation so it really can only beused in conjunction with the
latter boundary condition.

One caveat that must be mentioned is that the formulation of these two boundary
conditions is not general and therefore they should only be applied to web geometries
that are predominantly horizontal. That is, the x component of the normal vector to the
web sideset should at each point be less than or equal to the y component.



