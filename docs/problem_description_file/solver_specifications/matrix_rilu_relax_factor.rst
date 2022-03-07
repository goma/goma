****************************
Matrix RILU Relax Factor
****************************

::

	Matrix RILU relax factor = <float>

-----------------------
Description / Usage
-----------------------

This optional card provides a relaxation factor to Aztec to be used in conjunction with
preconditioners based on RILU(k,ω) approximate factorization. The input parameter
<float> is defined as

<float>
    **fac,** a floating point number ( fac ≥ 0 ) that specifies a relaxation
    factor.

If the *Matrix RILU relax factor* card is omitted, the default is **1.**

------------
Examples
------------

Following is a sample card:
::

	Matrix RILU relax factor = 0.5

-------------------------
Technical Discussion
-------------------------

Some limiting values for fac provide specific behavior:

* for a value of zero, the ILU(k) is obtained

* for a value of one, the MILU(k) is obtained.

The value of *k* is set by the *Matrix graph fillin* card.



