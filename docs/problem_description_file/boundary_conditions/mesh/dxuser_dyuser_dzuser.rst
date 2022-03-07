************************
**DXUSER DYUSER DZUSER**
************************

::

	BC = {DXUSER | DYUSER | DZUSER} SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(PCC/MESH)**

This boundary condition format is used to set a constant X, Y, or Z displacement as a
function of any independent variable available in Goma. These boundary conditions
require the user to edit the routines *dx_user_surf*, *dy_user_surf*, and/or *dz_user_surf* to
add the desired models. These routines are located in the file *user_bc.c.* In the input
deck each such specification is made on a separate input card. These boundary
conditions must be applied to side sets. Definitions of the input parameters are as
follows:

================================= =============================================================
**{DX_USER | DY_USER | DZ_USER}** Seven-character boundary condition
                                  name (<bc_name>) that defines the displacement, where:
                                 	
                                 	* **DX_USER** - X displacement, user-defined
                                 	* **DY_USER** - Y displacement, user-defined
                                 	* **DZ_USER** - Z displacement, user-defined
**SS**                            Type of boundary condition (<bc_type>), where **NS** denotes
                                  node set in the EXODUS II database.
<bc_id>                           The boundary flag identifier, an integer associated with
                                  <bc_type> that identifies the boundary location (side set
                                  in EXODUS II) in the problem domain.
<float_list>                      A list of float values separated by spaces which will be
                                  passed to the user-defined subroutine so the user can
                                  vary the parameters of the boundary condition. This list
                                  of float values is passed as a one-dimensional double
                                  array to the appropriate C function.
================================= =============================================================

------------
**Examples**
------------

Following is a sample card which applies an X-displacement boundary condition to the
nodes in node set 100, with a functional form set by the user and parameterized by the
single floating point number . These displacements are applied immediately to the
unknowns and hence result in immediate mesh displacement from the initial state.

::

     BC = DX_USER SS 100 1.0

Please consult the user-definition subroutines for examples. 




--------------
**References**
--------------

No References.