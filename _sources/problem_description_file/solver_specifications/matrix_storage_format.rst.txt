*************************
Matrix Storage Format
*************************

::

	Matrix storage format = {msr | vbr}

-----------------------
Description / Usage
-----------------------

This optional card can be used to choose between two formats accepted by the Aztec
2.1 solver package. Valid options are:

msr
    modified sparse row format (see Schunk and Shadid, 1992). This option is
    the default option and is automatically used for all direct solver options.
vbr
    variable block row format (see Heroux, 1992). This option should only be
    selected when an Aztec iterative solver is chosen.
epetra
    Compressed Sparse Row format using the Epetra library from Trilinos

------------
Examples
------------

Following is a sample card:
::

	Matrix storage format = msr

-------------------------
Technical Discussion
-------------------------

*Goma* supports two global matrix formats for its linear solvers. The advantage of
choosing **vbr** over the default **msr** format is a matter of which preconditioner option is
selected. (See Schunk, et al., 2002 on iterative methods.) When using the front solver
package, another format known as **estifm** is employed internally but not specified by
this card, which is not used in this case.

--------------
References
--------------

SAND92-1158: Iterative Solvers in Implicit Finite Element Codes, Sandia Technical
Report, Schunk, P. R. and Shadid, J. N. (1992)

SAND2001-3512J: Iterative Solvers and Preconditioners for Fully-coupled Finite
Element Formulations of Incompressible Fluid Mechanics and Related Transport
Problems, P. R. Schunk, M. A. Heroux, R. R. Rao, T. A. Baer, S. R. Subia and A. C.
Sun, March 2002.

TR/PA/92/90: M. A. Heroux, A proposal for a sparse BLAS toolkit, Technical Report,
CERFACS, December 1992.
