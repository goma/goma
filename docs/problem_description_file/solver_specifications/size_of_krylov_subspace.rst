***************************
Size of Krylov Subspace
***************************

::

	Size of Krylov subspace = <integer>

-----------------------
Description / Usage
-----------------------

This optional card allows the user to specify the dimension (size) of the Krylov
subspace for the **gmres** option of the *Solution Algorithm* card, where

<integer>
    **m,** specifies the number of orthogonalization directions and can be any
    positive integer less than or equal to the order of the matrix.

If the *Size of Krylov subspace* card is omitted, then the default dimension is **m** = 30.

------------
Examples
------------

The following is a sample input card:
::

	Size of Krylov subspace = 128

-------------------------
Technical Discussion
-------------------------

If the size of the subspace is at least as large as the maximum number of iterations
permitted by the solver then the **gmres** iteration will not include any restarts.
Depending on the problem, restarts may be beneficial, and then again they may not.
Particularly poorly conditioned linear systems may never converge below a certain
tolerance if **gmres** is allowed to restart (i.e. they “level off”). However, some linear
systems will admit a converged solution more rapidly with restarts than without.
Consequently, the user may wish to experiment with different values of this parameter.
See the Orthogonalization card for related information.

**gmres’** internal iterations create a Krylov subspace up to dimension m (less in some
circumstances, such as convergence). The time and space required by the internal
iterations increases nonlinearly with **m** (but see the *Orthogonalization* card) - a
doubling of **m** will result in more than a doubling of space and time requirements. So
simply choosing a very large dimension is generally not recommended.

