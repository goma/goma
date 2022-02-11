********************
Time Integration
********************

::

	Time integration = {steady | transient}

-----------------------
Description / Usage
-----------------------

This required card is used to specify transient or steady-state calculation. Valid options
are:

steady
    For a solution to the steady (time-derivative free) equations.

transient
    For transient simulations.

If option **steady** is chosen, then none of the other Time Integration Specification cards
in this section are needed.

------------
Examples
------------

This is a sample card for a steady state simulation:
::

	Time integration = steady

This is a sample card for a transient simulation:
::

	Time integration = transient

