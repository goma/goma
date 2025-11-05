***********************
Mesh correction damping
***********************

::

	Mesh correction damping = <float_list>

-----------------------
Description / Usage
-----------------------

<float1>
   Damping choices to be associated with Mesh correction tolerance.

------------
Examples
------------


Example :
::

	Newton line search type = BACKTRACK_MESH
    Mesh correction damping = 0.5 0.2 0.1
    Mesh correction tolerance = 1e-5 1e-4 1e-3

if :math:`L_2` norm of the mesh displacement residual > 1e-3, the damping factor is 0.1
if :math:`L_2` norm of the mesh displacement residual > 1e-4, the damping factor is 0.2
if :math:`L_2` norm of the mesh displacement residual > 1e-5, the damping factor is 0.5
else the damping factor is 1.0 (no damping)

This is applied only to mesh variables when using BACKTRACK_MESH option in Newton line search type card.

These are the current defaults if not specified.

-------------------------
Technical Discussion
-------------------------