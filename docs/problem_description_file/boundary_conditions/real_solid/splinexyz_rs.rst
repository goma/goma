****************
**SPLINEXYZ_RS**
****************

-----------------------
**Description / Usage**
-----------------------

**(PCC/MESH)**

This card is used to specify a general surface (solid) boundary description for
TOTAL_ALE real solid equations (see *Mesh Motion* card). These boundary conditions
are tantamount to *SPLINE_RS* , except that they do *not* invoke a vector residual
rotation into normal-tangential form. Instead, *SPLINEX_RS* invokes the geometric
boundary condition on the x-component of the real solid equation residual, and so on.
The card requires user-defined subroutines. Templates for these routines are currently
located in the routine “user_bc.c”. Both a function routine, fnc, for function
evaluation and corresponding routines dfncd1, dfncd2, and dfncd3 for the
derivative of the function with respect to global coordinates are required. Note that it
takes an arbitrary number of floating-point parameters, depending on the user’s needs.

Definitions of the input parameters are as follows:

+--------------+----------------------------------------------------------+
|{bc_name}     | Boundary condition name that defines the general surface;|
|              | the options are:                                         |
|              |                                                          |
|              |   * **SPLINEX_RS** - X general surface                   |
|              |   * **SPLINEY_RS** - Y general surface                   |
|              |   * **SPLINEZ_RS** - Z general surface                   |
+--------------+----------------------------------------------------------+
|**SS**        | Type of boundary condition (<bc_type>), where **SS**     |
|              | denotes side set in the EXODUS II database.              |
+--------------+----------------------------------------------------------+
|<bc_id>       | The boundary flag identifier, an integer associated with |
|              | <bc_type> that identifies the boundary location (side    |
|              | set in EXODUS II) in the problem domain.                 | 
+--------------+----------------------------------------------------------+
|[floatlist]   | Constants to parameterize any f(x,y,z) = 0 function      |
|              | input in user-defined routine fnc.                       |      
+--------------+----------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:
::

     BC = SPLINEZ_RS SS 10 1.0 100. 20.0 1001.0 32.0

applies a user-defined distinguishing condition parameterized by the list of floating
points to the boundary defined by side set 10. Most importantly, the condition replaces
the Z-component of the real solid equation.

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



--------------
**References**
--------------

No References.