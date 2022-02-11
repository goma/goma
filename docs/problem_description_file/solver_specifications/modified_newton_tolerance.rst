*****************************
Modified Newton Tolerance
*****************************

::

	Modified Newton Tolerance = <float1> <float2>

-----------------------
Description / Usage
-----------------------

This optional card allows the user to exert finer control over Jacobian formation than a
stride specification (as with the *Number of Newton Iterations* card’s second parameter
or the *Jacobian Reform Time Stride* card). Input parameters are defined as:

<float1>
    **r,** if the convergence rate is below this level ( r > 0.0 ), a Jacobian
    reformation will be forced.
<float2>
    **t,** if the residual norm is above this level ( t ≥ 0.0 ), a Jacobian
    reformation will be forced.

If the *Modified Newton Tolerance* card is omitted, then reformations are always
computed, subject to the *Number of Newton Iterations’* second parameter and the
*Jacobian Reform Time Stride* value.

See the *Jacobian Reform Time Stride* card for some detailed examples of the interaction
amongst various cards that influence when a Jacobian reformation occurs.

------------
Examples
------------

Following is a sample card:
::

	Modified Newton Tolerance = 1.5 1.0e-8

-------------------------
Technical Discussion
-------------------------

The convergence rate is defined as:

 .. math::

    \mathrm{convergence} \, \mathrm{rate} = \frac{\log \left( \mathrm{current} L_1 \mathrm{norm} \right) }{\log \left( \mathrm{previous} L_1 \mathrm{norm} \right)}

This rate should be equal to 2 when Newton’s method is in its region of convergence
(this is what it means to converge quadratically). A secant method would have a
convergence rate of :math:`1 + \sqrt{5}/2` (the golden ratio!), approximately 1.6.

The residual norm is simply the L\ :sub:`1` norm of the residual after a Newton iteration.

The method used to determine if a Jacobian reformation should take place is
conservative. If either test condition for reformation is satisfied, a reformation occurs.
Often, this card will allow you to speed up your runs by foregoing a fresh Jacobian reformation, but still maintain strong convergence. Moreover, without a Jacobian
reformation, the **lu** solver (see the *Solution Algorithm* card) can use a previously
factored matrix and simply do a resolve.
