***********
**UVWUSER**
***********

::

	BC = {UUSER | VUSER | WUSER} SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MOMENTUM)**

This card permits the user to specify an arbitrary integrated condition to replace a
component of the fluid momentum equations on a bounding surface. Specification of
the integrand is done via the functions uuser_surf, vuser_surf and
wuser_surf in file “user_bc.c.”, respectively.

A description of the syntax of this card follows:

+----------------------------+-------------------------------------------------------------+
|**{UUSER | VUSER | WUSER}** | Five-character boundary condition name (<bc_name>)          |
|                            | identifies the momentum equation component:                 |
|                            |                                                             |
|                            |   * **UUSER** - X momentum component                        |
|                            |   * **VUSER** - Y momentum component                        |
|                            |   * **WUSER** - Z momentum component                        |
+----------------------------+-------------------------------------------------------------+
|**SS**                      | Type of boundary condition (<bc_type>), where **SS**        |
|                            | denotes side set in the EXODUS II database.                 |
+----------------------------+-------------------------------------------------------------+
|<bc_id>                     | The boundary flag identifier, an integer associated with    |
|                            | <bc_type> that identifies the boundary location (side set   |
|                            | in EXODUS II) in the problem domain.                        |
+----------------------------+-------------------------------------------------------------+
|<float_list>                | A list of float values separated by spaces which will be    |
|                            | passed to the user-defined subroutines so the user can      |
|                            | vary the parameters of the boundary condition. This list    |
|                            | of float values is passed as a one-dimensional double       |
|                            | array to the appropriate C function.                        |
+----------------------------+-------------------------------------------------------------+


------------
**Examples**
------------

The following is an example of card syntax:
::

     BC = VUSER SS 10 1.0

Implementing the user-defined functions requires knowledge of basic data structures in
*Goma* and their appropriate use. The uninitiated will not be able to do this without
guidance.



--------------
**References**
--------------

No References.