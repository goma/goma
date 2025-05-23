***********************
Matrix Graph Fillin
***********************

::

	Matrix graph fillin = <integer>4.7.10_matrix_factorization_overlap.txt

-----------------------
Description / Usage
-----------------------

This optional card sets the graph level of fill-in for approximate factorizations used in
preconditioner construction for ILU(k), ICC(k) and BILU(k). The input parameter is
defined as

<integer>
    k, specifies the graph level of fill-in, k > 0.

If the *Matrix graph fillin* card is omitted, the default value of **k** is **0.**

------------
Examples
------------

Following is a sample card:
::

	Matrix graph fillin = 2

-------------------------
Technical Discussion
-------------------------

As the level of graph fill-in increases, the accuracy (usefulness) of the preconditioner
increases; however, so does memory usage as well as the time required to compute the
preconditioner.

