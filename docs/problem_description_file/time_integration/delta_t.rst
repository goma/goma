***********
delta_t
***********

::

	delta_t = <float>

-----------------------
Description / Usage
-----------------------

This card is required for **transient** simulations to set the value of the initial time step.
The input parameter is defined as:

<float>
    Any floating point number that indicates the time step in the appropriate
    units for your problem.

To specify a fixed time step size for an analysis, set <float> to be a negative number,
e.g. -1.0e-6; the code will use a constant (positive) time step. Should convergence
problems occur when a fixed step size is specified, the size of the time increment
entered for the *delta_t* card will be reduced by half until convergence is achieved. Once
a constant time step is reduced, it will not be increased.

------------
Examples
------------

Following is a sample card for an initial time step:
::

	delta_t = 6.e-03

If a constant time step is desired, use a negative value:
::

	delta_t = -6.e-03

