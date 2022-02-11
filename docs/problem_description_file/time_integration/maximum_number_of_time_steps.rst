********************************
Maximum Number of Time Steps
********************************

::

	Maximum number of time steps = <integer>

-----------------------
Description / Usage
-----------------------

This card sets the maximum number of time steps that may be performed for a
**transient** simulation. *Goma* will stop if this limit is reached. The input parameter is
defined as

<integer>
    Any integer greater than zero, which will limit the number of time steps
    taken in a simulation.

------------
Examples
------------

The following sample card sets the maximum number of time steps to 100:
::

	Maximum number of time steps = 100

