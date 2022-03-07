**********************
Matrix Output Type
**********************

::

	Matrix output type = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card indicates a level of diagnostic output for Aztec. The valid input
parameters for {char_string} are either a string or a positive integer:

all
    Print matrix and indexing vectors for each processor and all intermediate
    residual expressions.
none
    No intermediate results are printed. This is the default.
warnings
    Only Aztec warnings are printed.
last
    Only the final residual expression is printed.
k
    Residual expressions are printed every *k* iterations, k > 0.      

If the *Matrix output type* card is omitted, the default is **none.**

------------
Examples
------------

Following is a sample card:
::

	Matrix output type = 10

