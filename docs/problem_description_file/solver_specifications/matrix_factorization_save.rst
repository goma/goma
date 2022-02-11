*****************************
Matrix Factorization Save
*****************************

::

	Matrix factorization save = {0 | 1}

-----------------------
Description / Usage
-----------------------

This optional card is a boolean specification that determines whether the preconditioner
factorization information should be kept after a solve. Valid options are

0
    Factorization information is discarded.
1
    Factorization information is kept for that step.

If the *Matrix factorization save* card is omitted, then the default selection is **0.**

------------
Examples
------------

Following is a sample card:
::

	Matrix factorization save = 1

-------------------------
Technical Discussion
-------------------------

This option is most useful for iterative solution techniques where the computed
preconditioning matrix found from an incomplete factorization requires significant
computational resources. Such a preconditioner may be useful in later matrix solves
and obviate the need to compute another expensive preconditioner at the later stage.
Although a lot of time may be saved by re-using a previous factorization, the loss in
accuracy may cause convergence problems.



