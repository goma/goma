***********
**UVWVARY**
***********

::

	BC = {UVARY | VVARY | WVARY} SS <bc_id> [float_list]

-----------------------
**Description / Usage**
-----------------------

**(PCC/MOMENTUM)**

The *UVARY, VVARY* and *WVARY* boundary condition format is used to set variation in
X, Y, or Z velocity component, respectively, with respect to coordinates and time on a
specified sideset. Each such specification is made on a separate input card.

The *UVARY, VVARY*, and *WVARY* cards each require user-defined functions be supplied
in the file user_bc.c. Four separate C functions must be defined for a boundary
condition: velo_vary_fnc, dvelo_vary_fnc_d1, dvelo_vary_fnc_d2,
and dvelo_vary_fnc_d3. The first function returns the velocity component at a
specified coordinate and time value, the second, third, and fourth functions return the
derivative of the velocity component with x, y and z respectively.

A description of the syntax of this card follows:

+----------------------------+-------------------------------------------------------------+
|**{UVARY | VVARY | WVARY}** | Five-character boundary condition name (<bc_name>)          |
|                            | identifies the velocity component:                          |
|                            |                                                             |
|                            |   * **UVARY** -X velocity component                         |
|                            |   * **VVARY** -Y velocity component                         |
|                            |   * **WVARY** -Z velocity component                         |
+----------------------------+-------------------------------------------------------------+
|**SS**                      | Type of boundary condition (<bc_type>), where **SS**        |
|                            | denotes side set in the EXODUS II database.                 |
+----------------------------+-------------------------------------------------------------+
|<bc_id>                     | The boundary flag identifier, an integer associated with    |
|                            | <bc_type> that identifies the boundary location (side set   |
|                            | in EXODUS II) in the problem domain.                        |
+----------------------------+-------------------------------------------------------------+
|[float_list]                | An optional list of float values separated by spaces        |
|                            | which will be passed to the user-defined subroutines to     |
|                            | allow the user to vary the parameters of the boundary       |
|                            | condition. This list of float values is passed as a         |
|                            | onedimensional double array designated p in the parameter   |
|                            | list of all four C functions.                               |
+----------------------------+-------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card for an X component
::

     BC = UVARY SS 10 2.0 4.0

Following are the C functions that would have to be implemented in “user_bc.c” to
apply the preceding boundary condition card to set a parabolic velocity profile along a
sideset.

::

     double velo_vary_fnc( const int velo_condition, const double x,
       const double y, const double z, const double p[], const double
       time )
     {
     double f = 0;
     double height = p[0];
     double max_speed = p[1];
       if ( velo_condition == UVARY ) {
       f = max_speed*( 1.0 - pow(y/height, 2 ) );
       }
     return(f);
     }
     /* */
     double dvelo_vary_fnc_d1( const int velo_condition, const double
       x, const double y, const double z, const double p[], const
       double time )
     {
     double f = 0;
     return(f);
     }
     /* */
     double dvelo_vary_fnc_d2( const int velo_condition, const double
       x, const double y, const double z, const double p[], const
       double time )
     {
     double f = 0;
     double height = p[0];
     double max_speed = p[1];
       if ( velo_condition == UVARY ) {
       f = -2.0*max_speed*(y/height)/height;
     }
     return(f);
     }
     /* */
     double dvelo_vary_fnc_d3( const int velo_condition, const double
       x, const double y, const double z, const double p[], const
       double time )
     {
     double f = 0;
     return(f);
     }
     /* */

-------------------------
**Technical Discussion**
-------------------------

* Including the sensitivities is a pain, but required since *Goma* has no provision for
  computing Jacobian entries numerically.

* Note that the type of boundary condition (*UVARY, VVARY*, or *WVARY*) is sent to
  each function in the velo_condition parameter. Since there can be only one
  set of definition functions in user_bc.c, this allows the user to overload these
  functions to allow for more than one component defined in this manner. It would
  also be possible to use these functions to make multiple definitions of the same
  velocity component on different sidesets. However, this would have to be done by
  sending an identifier through the p array.

* This is a collocated-type boundary condition. It is applied exactly at nodal
  locations but has lower precedence of application than direct Dirichlet conditions.



--------------
**References**
--------------

No References.