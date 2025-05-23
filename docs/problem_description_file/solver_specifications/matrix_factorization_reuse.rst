******************************
Matrix Factorization Reuse
******************************

::

	Matrix factorization reuse = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card directs the approximate factorization solvers used in preconditioner
construction to reuse matrix information that may have been obtained during previous
linear solution stages. This card only has an effect when using an Aztec solver. Valid
options for {char_string} are:

calc
    Use no information from previous linear solutions.
recalc
    Use information from previous linear solutions but recalculate the
    preconditioning factors, with the implication that the symbolic
    factorization will be similar.
reuse
    Use information from previous linear solution; do not recalculate
    preconditioner factorizations. However, use scaling factors from previous
    linear solutions to scale righthand sides, initial guesses, and final
    solutions.

If the *Matrix factorization* reuse card is omitted, the default is **recalc.**

------------
Examples
------------

Following is a sample card:
::

	Matrix factorization reuse = recalc

-------------------------
**Technical Discussion**
-------------------------

See related discussions for *Matrix factorization save*.

