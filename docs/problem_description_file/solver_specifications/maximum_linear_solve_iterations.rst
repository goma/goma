***********************************
Maximum Linear Solve Iterations
***********************************

::

	Maximum Linear Solve Iterations = <integer>

-----------------------
Description / Usage
-----------------------

This optional card limits the maximum number of iterations used by iterative linear
solver algorithms. The input parameter is defined as

<integer>
    <integer> **n,** any positive integer ( n > 0 ) that specifies the maximum
    number of iterations.

If the *Maximum Linear Solve Iterations* card is omitted, the default selection is 500.

------------
Examples
------------

Following is a sample card:
::

	Maximum Linear Solve Iterations = 5

-------------------------
Technical Discussion
-------------------------

If the linear system can be solved within a specified tolerance (see the *Residual Ratio
Tolerance* card) in less than **n** iterations, then a normal return from Aztec occurs and
the actual number of iterations required to obtain convergence will be printed on the
status line. If the specified convergence tolerance is not met within **n** iterations, then an
abnormal return status occurs and, in place of the number of iterations, the string “max”
will be printed on the status line under the LIS (linear iteration status) heading. Other
abnormal returns from Aztec are possible and are indicated on the LIS status line; see
the Aztec User’s Guide (Hutchinson, Shadid and Tuminaro, 1995) for further
interpretation of different abnormal return status indicators.

--------------
References
--------------

SAND95-1559: Aztec User’s Guide Version 1.0, Sandia Internal Report, Hutchinson,
S. A., Shadid, J. N. and Tuminaro, R. S., 1995.
