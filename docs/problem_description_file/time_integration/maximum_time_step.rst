*********************
Maximum Time Step
*********************

::

	Maximum time step = <float>

-----------------------
Description / Usage
-----------------------

This card sets the value of the maximum allowable time step size in a **transient**
analysis, where the input parameter is defined as

<float>
    Any floating point number in the same units as specified in the *delta_t*
    card.

------------
Examples
------------

A sample card that sets the maximum time step to 10.0 follows:
::

	Maximum time step = 10.0

-------------------------
Technical Discussion
-------------------------

This setting is useful for advection dominated simulations, such as FILL, where a
Courant-like limit must be set on the value of the time step for optimal performance.

