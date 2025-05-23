******************************
Minimum Resolved Time Step
******************************

::

	Minimum Resolved Time Step = <float>

-----------------------
Description / Usage
-----------------------

Its role is to set a lower bound for the time step with respect to the **Time step error**
tolerance. When a converged time step is obtained by GOMA, the difference between
the predicted solution and final solution for that time step is compared to the **Time step
error** tolerance. If the difference exceeds this tolerance the step fails and the time step
is cut (usually by a factor of 2), UNLESS the time step falls below the **Minimum
Resolved Time Step** size. In this case the step is accepted, even if this error tolerance is
not achieved. This provides a mechanism for the modeler to control what phenomena is
resolved and what phenomena is ignored.

<float>
    Any floating point number in the same units as specified in the *delta_t*
    card.

------------
Examples
------------

A sample card that sets the maximum time step to 10.0 follows:
::

	Maximum Resolved Time Step = 10.0

-------------------------
Technical Discussion
-------------------------

See GT-034 for a thorough discussion.

--------------
References
--------------

GT-034: Tutorial on time step parameter selection for level-set problems in GOMA.
April 1, 2006. D. R. Noble
