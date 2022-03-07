***********************
Time Step Parameter
***********************

::

	Time step parameter = <float>

-----------------------
Description / Usage
-----------------------

This card allows the user to vary the time integration scheme. The usual settings are:

.. tabularcolumns:: |l|L|

=======================  ========================================================================
0.0                      Backward Euler method (1st order in time)
0.5                      Trapezoid rule (2nd order in time)
=======================  ========================================================================

------------
Examples
------------

This is a sample card that sets the time integration scheme to Trapezoidal rule:
::

	Time step parameter = 0.5

-------------------------
Technical Discussion
-------------------------

One should usually use the Trapezoid rule. When a large time step :math:`\Delta t` is used the
Trapezoid rule can exhibit oscillations. If such a large :math:`\Delta t` is required then the Backward
Euler method can be used (it will damp oscillations), albeit at a cost of accuracy.

If we designate the time step parameter as :math:`\theta`, the solution at time step :math:`n` as :math:`y^n`, and the
PDE to be solved as

.. math::

   \frac{\partial y}{\partial t} = g \left( y \right)

then the time integration method takes the form

.. math::

   \frac{y^{n+1} - y^n}{\Delta t} = \frac{2 \theta}{1 + 2 \theta} \dot{y}^n + \frac{1}{1 + 2 \theta} g \left( y^{n+1} \right)

where

.. math::

   \dot{y}^{n+1} = \frac{1 + 2 \theta}{\Delta t} \left( y^{n+1} - y^n \right) - 2 \theta \dot{y}^n = g \left( y^{n+1} \right).


Note that there is no choice of finite :math:`\theta` that will yield a Forward Euler method. See Gartling (1987) for more information.


--------
FAQs
--------

For porous flow problems with mass lumping, you should always choose backward
Euler method.

--------------
References
--------------

SAND86-1816: NACHOS 2: A Finite Element Computer Program for Incompressible
Flow Problems - Part 2 - User’s Manual, Gartling, David K. (September, 1987).

