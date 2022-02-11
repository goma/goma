***************************
Matrix Auxiliary Vector
***************************

::

	Matrix auxiliary vector = {resid | rand}

-----------------------
Description / Usage
-----------------------

This optional card indicates to Aztec how the auxiliary vector r is determined.
Permissible options are:

resid
    The auxiliary vector is set to the initial residual vector, viz. :math:`r = r(0)`.
rand
    The auxiliary vector is filled with random numbers, each in the range
    :math:`[-1,1]`.

If the *Matrix auxiliary vector* card is omitted, the default is **resid.**

------------
Examples
------------

Following is a sample card:
::

	Matrix auxiliary vector = rand

-------------------------
Technical Discussion
-------------------------

The auxiliary vector is only used for certain iterative linear matrix solution algorithms.

The **rand** option may cause difficulties with initial iterative solver steps because
different processors may have different initial unknown values at shared unknowns.
