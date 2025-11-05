***********************
Newton line search type
***********************

::

	Newton line search type = {FULL_STEP | BACKTRACK | BACKTRACK_MESH}

-----------------------
Description / Usage
-----------------------

This optional card has a single input parameter:

FULL_STEP
    A full step Newton iteration will be taken, this still can be controlled by manual relaxation such as through Newton correction factor

BACKTRACK
    Backtracking line search is used. The relaxation parameter is automatically chosen.

BACKTRACK_MESH
    Backtracking line search is used. The relaxation parameter is automatically chosen for all non-mesh variables. 
    The mesh variables are relaxed using a mesh correction damping and tolerance.

Default: Newton line search type = FULL_STEP

------------
Examples
------------


Example :
::

	Newton line search type = BACKTRACK

Example :
::

	Newton line search type = BACKTRACK_MESH
    Mesh correction damping = 0.5 0.2 0.1
    Mesh correction tolerance = 1e-5 1e-4 1e-3


-------------------------
Technical Discussion
-------------------------

Testing is still in progress and is subject to change. 

We use a simple algorithm for backtracking line search. The algorithm is as follows:

1. Compute the Newton step direction, p, by solving J * p = -R
2. Set the initial step length, alpha = 1
3. While the residual norm ||R(u + alpha * p)|| is not sufficiently decreased:
    a. Reduce the step length, alpha = alpha * 0.5
4. Update the solution, u = u + alpha * p
5. Repeat the step 3-4 until minimum step is reached or sufficient decrease is achieved.

