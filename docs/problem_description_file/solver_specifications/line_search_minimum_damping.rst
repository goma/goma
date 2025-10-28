***************************
Line search minimum damping
***************************

::

	Line search minimum damping = <float>

-----------------------
Description / Usage
-----------------------

<float>
   The minimum damping factor to be used in the backtracking line search
   algorithm. This is a value between 0 and 1.0. A value of 1.0 means no
   damping, while a value closer to 0 means more damping. This parameter is only
   used when the Newton line search type is set to BACKTRACK or BACKTRACK_MESH.

------------
Examples
------------


Example :
::

	Newton line search type = BACKTRACK_MESH
    Line search minimum damping = 0.01
    Mesh correction damping = 0.5 0.2 0.1
    Mesh correction tolerance = 1e-5 1e-4 1e-3
