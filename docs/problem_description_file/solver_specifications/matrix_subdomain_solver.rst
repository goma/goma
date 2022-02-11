***************************
Matrix Subdomain Solver
***************************

::

	Matrix subdomain solver = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card selects a solver to use in constructing a preconditioner. It is used in
conjunction with a *Preconditioner* setting.

::

	Preconditioner = dom_decomp

All of these preconditioners are available through the Aztec library. Valid options for
{char_string} are listed below.

lu
    Approximately solve the processor’s local matrix via direct factorization
    using Sparse 1.3 in conjunction with a userspecified *Matrix drop
    tolerance.*
ilut
    Approximately solve the processor’s local matrix via ILUT (Saad, 1994.) The
    factorization is affected by user-specified options for *Matrix drop
    tolerance* as well as *Matrix ILUT fill factor.*

    This subdomain solver is among the more robust to recommend as a first
    attempt; thus it has been chosen as the default if no subdomain solver is
    specified.
ilu
    Approximately factor the processor’s local matrix using ILU(k), where k is
    specified by the user in the argument to *Matrix graph fillin.*
rilu
    Approximately factor the processor’s local matrix using RILU(k,ω), where
    k is specified by the user in the argument to *Matrix graph fillin* and
    ω is specified by the user in the argument to *Matrix RILU relax factor.*
    (This option applies only to Trilinos.)
bilu
    Approximately factor the processor’s local matrix using block ILU(k) for
    a VBR format matrix, where k is specified by the user in the argument to
    *Matrix graph fillin.* While not the most efficient preconditioner,
    **bilu** is very robust. (This option applies only to Trilinos.)
icc
    Incomplete Cholesky factorization. See the Aztec manual for a reference.

If this *Matrix subdomain solver* card is omitted, then the default selection is **ilut.**

------------
Examples
------------

Following is a sample card:
::

	Matrix subdomain solver = ilut

-------------------------
Technical Discussion
-------------------------

There is no real recipe to follow when choosing a preconditioner. In general, the
cheapest preconditioner that works should be used. If ILUT(1) does the job, great.
Sometimes the only preconditioner(s) that will work are very expensive. When the
preconditioner seems to take too much time, remember that you may not be choosing
the “wrong” preconditioner; the problem may just be that difficult.

Although Aztec 2.1 is being maintained and supported as a solver package for Goma,
the interative solvers and preconditioners are now primarily accessed through the
Trilinos library (as AztecOO), which is actively being developed and maintained at
Sandia National Laboratories. Note that some features can only be accessed through
Trilinos, as indicated above.

--------------
References
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.

Saad, Y., 1994. “ILUT: a dual threshold incomplete ILU factorization”, Numerical
Linear Algebra with Applications, 1:387-402.
