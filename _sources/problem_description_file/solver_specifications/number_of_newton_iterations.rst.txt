*******************************
Number of Newton Iterations
*******************************

::

	Number of Newton Iterations = <integer1> [integer2]

-----------------------
Description / Usage
-----------------------

This required card sets the maximum number of iterations allowed for convergence of
the Newton nonlinear iteration loop. It also provides an optional parameter for setting
the reformation stride for the Jacobian matrix. Definitions of the input parameters are
as follows:

<integer1>
    n\ :sub:`1`, any integer indicating the maximum number of iterations
    allowed for convergence of the Newton nonlinear iteration loop.
[integer2]
    n\ :sub:`2`, an optional parameter indicating the reformation stride for
    the Jacobian matrix.

The *Number of Newton Iterations* card is required, there is no default.

See the *Jacobian Reform Time Stride* card for some detailed examples of the interaction
amongst various input parameters that influence when a Jacobian reformation occurs.

------------
Examples
------------

Following is a sample card:
::

	Number of Newton Iterations = 5

-------------------------
Technical Discussion
-------------------------

For an unrelaxed Newton iteration with a good initial guess, five or six iterations (for
n\ :sub:`1`) should be sufficient to achieve convergence for most problems. One iteration 
will
suffice for problems that are linear; two can be specified, with the second iteration
verifying that the residual norms are small. More iterations may be required for relaxed
Newton iteration schemes using the correction factor described in the *Newton
correction factor* card. This parameter can also be controlled from the command line
(see the **-n** option in the section on Command-line Arguments, Chapter 3).

The optional second parameter can be used to invoke a modified Newton iteration. If
this value is missing, the stride is set to unity. This capability enables the user to save
on assembly time when near a solution, particularly when doing transient simulations.

