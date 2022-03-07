****************
Maximum Time
****************

::

	Maximum time = <float>

-----------------------
Description / Usage
-----------------------

This card sets the maximum value of time that may be achieved in a **transient**
simulation. *Goma* will stop if this limit is reached. The input parameter is defined as:

<float>
    Any floating point number in the same units as specified in the *delta_t*
    card.

The last result written to the EXODUS II and soln.dat file in a successfully
completed simulation will always be at the maximum time. This provides a cutoff time
beyond which the simulation will terminate.

------------
Examples
------------

The following sample card sets the maximum time to 105 (in units consistent with your
simulation):
::

	Maximum time = 105.

