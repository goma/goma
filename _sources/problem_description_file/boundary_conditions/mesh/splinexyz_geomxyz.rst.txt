*********************
**SPLINEXYZ/GEOMXYZ**
*********************

::

	BC = {bc_name} SS <bc_id> [floatlist]

-----------------------
**Description / Usage**
-----------------------

**(PCC/MESH)**

This card is used to specify a general surface (solid) boundary description for ALE (or
in special cases LAGRANGIAN) type mesh motion (see *Mesh Motion* card). These
boundary conditions are tantamount to *SPLINE* or *GEOM*, except that they do *not*
invoke a mesh-equation vector residual rotation into normal-tangential form. Instead,
*SPLINEX* or, equivalently, *GEOMX* invokes the geometric boundary condition on the
x-component of the mesh equation residual, and so on. The card requires user-defined
subroutines. Templates for these routines are currently located in the routine
“user_bc.c”. Both a function routine, fnc, for function evaluation and
corresponding routines dfncd1, dfncd2, and dfncd3 for the derivative of the
function with respect to global coordinates are required. *GEOMX* and *SPLINEX* are
exactly the same condition. SPLINE* usage is being deprecated. Note that it takes an
arbitrary number of floating-point parameters, depending on the user’s needs.

Definitions of the input parameters are as follows:

================= ============================================================================
{bc_name}         Boundary condition name that defines the general surface;
                  the options are:

                  	* **SPLINEX/GEOMX** - X general surface
                  	* **SPLINEY/GEOMY** - Y general surface
                  	* **SPLINEZ/GEOMZ** - Z general surface
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

The following is a sample input card:
::

     BC = GEOMZ SS 10 1.0 100. 20.0 1001.0 32.0

applies a user-defined distinguishing condition parameterized by the list of floating
points to the boundary defined by side set 10. Most importantly, the condition replaces
the Z-component of the momentum equation.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form of this distinguishing condition is arbitrary and is specified by
the user in the fnc routine in user_bc.c. Derivatives of the user-specified function
must also be provided so as to maintain strong convergence in the Newton iteration
process. These functions are located next to fnc and are named dfncd1, dfncd2, and
dfncd3.Several examples for simple surfaces exist in the template routine. In three
dimensions, usage needs to be completed with a companion *ROT* input card which
directs the equation application of the condition, even though rotations are not actually
performed.



