**********************
Printing Frequency
**********************

::

	Printing Frequency = <integer> [float]

-----------------------
Description / Usage
-----------------------

This card sets the printing frequency, the step or time interval, at which *Goma* will print
the solution variables to the *Output EXODUS II file* and the *SOLN file*. Definitions of
the <integer> options, and the dependent [float] option when <integer> is set to 0, are:

<integer>
    Specifies how often the solution will be printed.

    .. tabularcolumns:: |l|L|

    ====   =====
    > 0    Interval in time steps between successive printings of the solution,
           any positive integer value.
    0      Controls printing of the solution at regularly spaced (uniform) intervals
           of time (every [float]), regardless of the number of time steps over that
           time interval
    ====   =====

[float]
    Elapsed time (in the same units as specified in the *delta_t* card) between
    successive printings of the solution (any positive number).

------------
Examples
------------

*Goma* will print the solution every five time steps given the following sample card:
::

	Printing Frequency = 5

*Goma* will print the solution every ten time units given the following sample card:
::

	Printing Frequency = 0 10.

