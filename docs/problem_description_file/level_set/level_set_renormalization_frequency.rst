***************************************
Level Set Renormalization Frequency
***************************************

::

	Level Set Renormalization Frequency = <integer>

-----------------------
Description / Usage
-----------------------

This card sets an upper limit to the number of time steps which are allowed to pass
between renormalization procedures. Possible values for <integer> are listed below:

<integer>

    .. tabularcolumns:: |l|L|

    ==== ====
    -1   never renormalize (default)
    0    renormalize every step
    *n*  a positive integer >1, renormalize every n\ :sup:`th` time step
    ==== ====

------------
Examples
------------

This is a sample input:
::

	Level Set Renormalization Frequency = 50

-------------------------
Technical Discussion
-------------------------

Renormalization procedures are normally triggered by the average gradient exceeding
one by a specified amount (*see Level Set Renormalization Tolerance*). However, at
times it might be advantageous to trigger a renormalization independent of the size of
the average level set gradient. For example, it might occur that in a very small region
near the interface, the level set gradient is becoming large but elsewhere the gradient is
still relatively small. Since the average gradient is used, this condition might not trigger
renormalization. By setting an upper limit for the number of time steps that can pass
before renormalization, situations such as this can be remedied.



