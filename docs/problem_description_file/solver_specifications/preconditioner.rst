******************
Preconditioner
******************

::

	Preconditioner = {char_string}

-----------------------
Description / Usage
-----------------------

Iterative techniques for solving a linear matrix system (see above) often benefit from
preconditioning to aid convergence. This optional card provides for the selection of a
preconditioner from those available through Aztec. For direct factorization *Solution
Algorithm* specifications, the *Preconditioner* specification is immaterial since none is
performed; in such cases, this card should be omitted.

Valid options for {char_string} are listed below.

none
    No preconditioning is performed. This is the default specification if no
    preconditioner has been specified.
Jacobi
    A k-step Jacobi preconditioner is used (block Jacobi for VBR matrices). The
    number of Jacobi steps, k, is set using the *Matrix polynomial order* card.
Neumann
    A Neumann series polynomial preconditioner is used, where the order of the
    polynomial, k, is set using the *Matrix polynomial order* card.
ls
    A least-squares polynomial preconditioner is used, where the order of the
    polynomial, k, is set using the *Matrix polynomial order* card.
sym_GS
    A k-step symmetric Gauss-Seidel preconditioner is used for non-overlapping
    domain decomposition (additive Schwarz). In parallel, each processor
    performs one step of symmetric Gauss-Seidel on its local matrix, followed
    by communication to update boundary values from adjacent processors before
    performing the next local symmetric Gauss-Seidel step. The number of steps,
    k, is set using the *Matrix polynomial order* card.
lu
    Approximately solve the processor’s local matrix via direct factorization
    using Sparse 1.3 in conjunction with a userspecified *Matrix drop
    tolerance.*
dom_decomp
    A domain-decomposition-based preconditioner (additive Schwarz). Each
    processor augments its local matrix according to the *Matrix factorization
    overlap* card and then approximately solves the resulting linear system
    using the solver specified by the *Matrix subdomain solver* card. This is
    the most often used *Preconditioner* card.

------------
Examples
------------

Following is a sample card:
::

	Preconditioner = dom_decomp

-------------------------
Technical Discussion
-------------------------

Note that prior to Aztec 2.x, certain subdomain solvers were specified simply as
arguments to the *Preconditioner* card. While this historical usage is permitted via
limited backward compatibility in order to ease the transition from Aztec 1 usage, the
preferred usage is to specify ILU (and similar) preconditioners as a subdomain solver
using the more powerful and flexible options that are available using Aztec 2.x together
with this option for the preconditioner. Since subdomain solvers such as ILU and ILUT
are powerful and frequently used, this preconditioner option will predominate when
iterative solvers are being used, even in serial execution.

The most popular setting is **dom_decomp**, with a subdomain solver specified in the
*Matrix Subdomain Solver* card. For further details, consult Mike Heroux’s recipe for
applying preconditioners and what to dial the knobs to (in Schunk, et. al., 2002).

--------------
References
--------------

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.
