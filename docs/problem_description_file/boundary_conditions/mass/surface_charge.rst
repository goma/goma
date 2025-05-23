******************
**SURFACE_CHARGE**
******************

::

	BC = SURFACE_CHARGE SS <bc_id> <integer> <float>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MASS)**

The *SURFACE_CHARGE* card specifies the electrostatic nature of a surface:
electrically neutral, positively charged or negatively charged.

Definitions of the input parameters are as follows:

+------------------+-----------------------------------------------------------+
|**SURFACE_CHARGE**| Name of the boundary condition (<bc_name>).               |
+------------------+-----------------------------------------------------------+
|**SS**            | Type of boundary condition (<bc_type>), where **SS**      |
|                  | denotes side set in the EXODUS II database.               |
+------------------+-----------------------------------------------------------+
|<bc_id>           | The boundary flag identifier, an integer associated with  |
|                  | <bc_type> that identifies the boundary location (side set |
|                  | in EXODUS II) in the problem domain.                      |
+------------------+-----------------------------------------------------------+
|<integer>         | Index of species to which surface charge condition        |
|                  | applies.                                                  |
+------------------+-----------------------------------------------------------+
|<float>           | *z*, value of surface charge.                             |
+------------------+--------------+--------------------------------------------+
|                  | 0            | electroneutrality                          |
+------------------+--------------+--------------------------------------------+
|                  | positive *z* | positively charged surface                 |
+------------------+--------------+--------------------------------------------+
|                  | negative *z* | negatively charged surface                 |
+------------------+--------------+--------------------------------------------+

------------
**Examples**
------------

The following input card indicates that on side set 1 species 1 is electrically neutral:
::

   BC = SURFACE_CHARGE SS 1   1   0.0
