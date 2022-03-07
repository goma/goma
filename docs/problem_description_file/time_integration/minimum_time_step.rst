*********************
Minimum Time Step
*********************

::

	Minimum time step = <float>

-----------------------
Description / Usage
-----------------------

This card sets the value of the minimum allowable time step size in a **transient**
analysis, a useful control if the time step is being decreased due to poor convergence of
the transient or iterative algorithm. The input parameter is defined as

<float>
    Any floating point number in the same units as specified in the *delta_t*
    card.

------------
Examples
------------

A sample card that sets the minimum time step to 1.e-9 follows:
::

	Minimum time step = 1.e-9

-------------------------
Technical Discussion
-------------------------

This specification provides a graceful way for the program to terminate based on the
computed time step dropping below the minimum value rather than terminating by a
segmentation fault or a divide-by-zero error that could result if the time step becomes
too small without the benefit of this control.

