*******************
Time Step Error
*******************

::

	Time step error = <float> <integer_list>

-----------------------
Description / Usage
-----------------------

The time step error controls the adjustable time step size based on the difference
between the solution and the predicted solution (L\ :sub:`2` norm). The first of the eight
arguments is a floating point number that indicates the error in the time step selection.

<float>
    the error value, any floating point number.

    The smaller this number is, the smaller the time step will tend to be in the automatic
    time step control. The original implementation of this capability in *Goma* did not use a
    normalized value for the norm; to enable this most useful feature, use a negative value
    of the time step error and a positive, normalized norm will be computed. This way a
    percentage value of the solution error will be set.

<integer_list>
    seven integers, with a value either zero (0) or one (1).

    A further degree of control is offered by the seven integers (i\ :sub:`1`
    through i\ :sub:`7`) that identify which solution variables will contribute
    to the error norm calculations. Permissible values for each of these seven
    integers are 0 and 1. The correspondence between the integers and variables
    is as follows:

    .. tabularcolumns:: |l|L|

    =======================  ========================================================================
    i\ :sub:`1`              (pseudo) solid displacement
    i\ :sub:`2`              fluid velocity
    i\ :sub:`3`              temperature
    i\ :sub:`4`              concentration, porous liquid pressure, gas pressure, porosity,
                             saturation
    i\ :sub:`5`              pressure
    i\ :sub:`6`              fluid (polymer) extra stress
    i\ :sub:`7`              voltage
    =======================  ========================================================================

    A value of 0 for an integer directs *Goma* to exclude contributions from
    that variable in the error norm calculation; correspondingly, a value of
    **1** means that variable should be included.

------------
Examples
------------

A sample time step error card follows:
::

	Time step error = 0.01 0 1 1 1 0 0 0

In this example, the L\ :sub:`2` norms for the fluid velocity, temperature, and concentration are
summed (and scaled) prior to comparison with the target error value of 0.01. If the
norms of the velocity, temperature, and concentration variables is greater than 0.01, the
time step is halved and the step repeated. Otherwise, the current step size is compared
to other step criteria before continuing to the next step.

If the integer values are omitted, the scaled error norm becomes infinite and the
analysis will terminate in the error norm calculation with an arithmetic overflow.

------------
Examples
------------

To use the normalized value of the norm, the following would be specified:
::

	Time step error = -0.01 0 1 1 1 0 0 0

This would set the maximum time step error to be 1%.

-------------------------
Technical Discussion
-------------------------

Note that on porous flow problems the error in step-size is computed as a composite
measure of all porous-flow variables, viz. these cannot currently be controlled
separately.

