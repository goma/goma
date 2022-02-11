*******************************
Jacobian Reform Time Stride
*******************************

::

	Jacobian Reform Time Stride = <integer>

-----------------------
Description / Usage
-----------------------

This optional card has a single input parameter:

<integer>
    **k,** the stride length for Jacobian reformations ( :math:`k \geq 1` ).

The *Jacobian Reform Time Stride* card is optional; there is no default.

------------
Examples
------------

Three examples are provided to illustrate how to use this card.

Example #1:
::

	Number of Newton Iterations = 12 1

	Modified Newton Tolerance = 1.9 0.1

	Jacobian Reform Time Stride = 2

	Newton correction factor = 1

This will reform the Jacobian every 2 steps. Furthermore, if the convergence rate falls
below 1.9 or the L\ :sub:`1` residual is greater than 0.1 on an off-stride step a Jacobian reformation
will occur. Specifically, the*Modified Newton Tolerance* takes precedence over a
reformation stride setting (from either *Number of Newton Iterations or Jacobian
Reform Time Stride*).

Example #2:
::

	Number of Newton Iterations = 12 1

	# Modified Newton Tolerance = 1.9 0.1

	Jacobian Reform Time Stride = 2

	Newton correction factor = 1

Note this differs from the previous example only by omitting the *Modified Newton
Tolerance* card. This causes the Jacobian to be reformed every other time step.

Example #3:
::

	Number of Newton Iterations = 12 2

	# Modified Newton Tolerance = 1.9 0.1

	Jacobian Reform Time Stride = 1

	Newton correction factor = 1

We’ve changed the *Jacobian Reform Time Stride* from 2 to 1 and changed the second
parameter of the *Number of Newton Iterations* card from 1 to 2. This will cause the
Jacobian to be reformed every other step.

-------------------------
Technical Discussion
-------------------------

If the second parameter on the *Number of Newton Iterations* card is present and greater
than 1, this *Jacobian Reform Time Stride* card is ignored. Otherwise, this card simply
forces the Jacobian to be rebuilt every **k** Newton steps. Often, this card will allow you
to speed up your runs by foregoing a fresh Jacobian formation, but still maintain strong
convergence. Moreover, without a Jacobian formation, the **lu** solver (see the *Solution
Algorithm* card) can use a previously factored matrix and simply do a resolve.

