**********************************
Pressure Stabilization Scaling
**********************************

::

	Pressure Stabilization Scaling = <float>

-----------------------
Description / Usage
-----------------------

This optional card is only used if the Pressure Stabilization card is set to **yes,** where

<float>
    **tau,** a positive real value ( tau > 0.0 ) that scales the momentum
    residual being added to the continuity equation for pressure stabilization.

The default value of **tau** is 0.1. If the *Pressure Stabilization* card is omitted, or set to
**no,** then **tau** is ignored.

------------
Examples
------------

Following is a sample card:
::

	Pressure Stabilization Scaling = 0.01

-------------------------
Technical Discussion
-------------------------

Generally, if **tau** is small, then more accurate solutions may be obtained at the cost of a
more ill-conditioned matrix system that may not be easily amenable to iterative solvers
(but stay tuned!). Conversely, larger values of this parameter result in equation systems
that are easier to solve using the available iterative matrix solvers, but the solution thus
obtained may be less accurate. A good choice for **tau** is 0.1.

The scaling value, **tau,** is further scaled inside of *Goma*. Knowledge of this scaling is
sometimes useful. First, an average Reynolds number (*Re*) is computed according to:

.. math::

   \mathrm{Re} = \frac{\rho \lVert U \rVert \langle h \rangle }{2 \mu}

where :math:`\rho` and :math:`\mu` are local values for density and viscosity, :math:`\lVert U \rVert` is a norm of the velocity
field, and :math:`\langle h \rangle` is a global average value for element size. If Re < 3.0, the pressure
stabilization scaling is given by this expression:

.. math::

   \frac{ \mathrm{tau} {\langle h \rangle}^2 }{12 \mu}

On the other hand, if Re > 3.0, the following scales the pressure stabilization terms in
the continuity equation:

.. math::

   \frac{ \mathrm{tau} {\langle h \rangle} }{2 \rho \lVert U \rVert}   
