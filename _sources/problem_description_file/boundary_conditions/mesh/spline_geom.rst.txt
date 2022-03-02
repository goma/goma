***************
**SPLINE/GEOM**
***************

::

	BC = {SPLINE|GEOM} SS <bc_id> [floatlist]

-----------------------
**Description / Usage**
-----------------------

**(PCC/ROTATED MESH)**

This card is used to specify a general surface (solid) boundary description for ALE (or
in special cases LAGRANGIAN) type mesh motion (see *Mesh Motion* card). Like most
other distinguishing conditions, this condition causes the mesh-motion equations, viz.
*mesh1, mesh2,* and *mesh3*, to be rotated into boundary normal-tangential form. The
card requires user-defined subroutines. Templates for these routines are currently
located in the routine “user_bc.c”. Both a function routine, fnc, for function
evaluation and corresponding routines dfncd1, dfncd2, and dfncd3 for the
derivative of the function with respect to global coordinates are required. The SPLINE
condition is exactly the same and uses the same routine as the GEOM card option, and
hence as of the time of this writing we are deprecating the use of SPLINE. Note that it
takes an arbitrary number of floating-point parameters, depending on the user’s needs.

Definitions of the input parameters are as follows:

================= ============================================================================
**SPLINE/GEOM**   Name of the boundary condition <bc_name>).
**SS**            Type of boundary condition (<bc_type>), where **SS** denotes
                  side set in the EXODUS II database.
<bc_id>           The boundary flag identifier, an integer associated with
                  <bc_type> that identifies the boundary location (side set in
                  EXODUS II) in the problem domain.
[floatlist]       Constants to parameterize any f(x,y,z) = 0 function input in
                  user-defined routine fnc.
================= ============================================================================

------------
**Examples**
------------

The following sample input card:
::

     BC = SPLINE SS 10 1.0 100. 20.0 1001.0 32.0

applies a user-defined distinguishing condition, parameterized by the list of five
floating point values, to the boundary defined by side set 10.

-------------------------
**Technical Discussion**
-------------------------

This condition, like *DISTNG, PLANE*, and others that can be applied to geometry, is
applied to the normal component of the mesh motion equations along a boundary in
two dimensions; in three dimensions application needs to be further directed with the
*ROT* conditions. Examples of typical distinguishing conditions can be found in
user_bc.c in the fnc routine and companion derivative routines.




