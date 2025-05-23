******************
Matrix Scaling
******************

::

	Matrix scaling = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card selects a scaling for the linear matrix system solution step.Valid
options for {char_string} are listed below.

none
    No scaling is performed. This is the default if no *Matrix Scaling* card is
    present.
Jacobi
    Point Jacobi scaling is performed.
BJacobi
    Block Jacobi scaling is performed if the underlying matrix format is VBR.
    If the MSR matrix format is used, the scaling reverts to point Jacobi.
row_sum
    Scale each row so the sum of the magnitudes of the nonzero elements is 1.
sym_diag
    Symmetric scaling so that diagonal elements are 1.
sym_row_sum
    Symmetric scaling using the matrix row sums.

If the *Matrix Scaling* card is omitted, the default selection is **none.**

------------
Examples
------------

Following is a sample card:
::

	Matrix scaling = sym_diag

-------------------------
Technical Discussion
-------------------------

All of these scalings are supplied via the Aztec library and thus will not affect the linear
systems that are solved by other means (using **front**, for example). In an odd twist of
fate, the linear system always undergoes a row sum scaling (equivalent to the **row_sum**
option) before these other scalings are applied. Note that when a nontrivial scaling is
selected, the matrix is overwritten with a rescaled system.

