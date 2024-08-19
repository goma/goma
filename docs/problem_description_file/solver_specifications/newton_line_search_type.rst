***********************
Newton line search type
***********************

::

	Newton line search type = {FULL_STEP | BACKTRACK}

-----------------------
Description / Usage
-----------------------

This optional card has a single input parameter:

FULL_STEP
    A full step Newton iteration will be taken, this still can be controlled by manual relaxation such as through Newton correction factor

BACKTRACK
    Backtracking line search is used. The relaxation parameter is automatically chosen.


Default: Newton line search type = FULL_STEP

------------
Examples
------------


Example :
::

	Newton line search type = BACKTRACK


-------------------------
Technical Discussion
-------------------------

Testing is still in progress and is subject to change. Currently backtracking line search does not work well for ALE problems that are likely to need manual relaxation.
