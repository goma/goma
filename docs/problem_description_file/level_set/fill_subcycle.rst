*****************
Fill Subcycle
*****************

::

	Fill Subcycle = <integer>

-----------------------
Description / Usage
-----------------------

This is an optional card that sets the number of subcycle-fill time steps between fluidflow
time steps in uncoupled level set calculations. The default is 10 subcycle time
steps for every flow time step. The input parameter is defined as

<integer>
    Any nonzero number indicating the subcycling frequency of the fill equation
    versus the flow equations.

For example, if the value of <integer> is 1, the flow and fill equations are solved every
time step. If it is 10, between every transient step in the flow calculation, the fill
(advection) equation is solved 10 times with one-tenth of the time step.

------------
Examples
------------

The following is a sample card that sets the fill subcycling rate to 4:
::

	Fill Subcycle = 4

--------------
References
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, Februray 27, 2001, T.A.
Baer
