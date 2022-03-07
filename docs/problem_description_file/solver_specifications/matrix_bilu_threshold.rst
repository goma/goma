*************************
Matrix BILU Threshold
*************************

::

	Matrix BILU threshold = <float>

-----------------------
Description / Usage
-----------------------

This capability is only present within the Trilinos library. This optional card provides a
means to modify the way the block ILU preconditioner (*Matrix subdomain solver* =
**bilu**) is constructed. The input parameter is defined as:

<float>
    **t,** a floating point number ( t ≥ 0.0 ) that sets the *Matrix Relative
    Threshold* and *Matrix Absolute Threshold* thresholds.

When the *Matrix BILU threshold* card is omitted, the default value is 0.0.

------------
Examples
------------

Following is a sample card:
::

	Matrix BILU Threshold = 1.0e-14

-------------------------
Technical Discussion
-------------------------

Using this card is equivalent to supplying both the *Matrix Relative Threshold* and
*Matrix Absolute Threshold* with the value specified with this card.

The value of **t** defaults to zero, and if given a small value, say 1.0e-14, the condition
number of the preconditioner, as reported when using the **bilu** option, should decrease.
Try increasing up to around 1.0e-3 to get added benefit. The **bilu** preconditioner is not
actually the cheapest or most efficient preconditioner, but it is very robust.



