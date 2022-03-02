***************************
Matrix ILUT Fill Factor
***************************

::

	Matrix ILUT fill factor = <float>

-----------------------
Description / Usage
-----------------------

This optional card provides a second criterion to Aztec to be used in conjunction with
preconditioners based on ILUT approximate factorization, where

<float>
    **fac,** a floating point value ( fac ≥ 0 ) that specifies, very crudely,
    how many nonzero entries the approximate factorization will contain
    relative to the number of nonzero entries in the original matrix.

If the *Matrix ILUT fill factor* card is omitted, the default is **1.**

------------
Examples
------------

Following is a sample card:
::

	Matrix ILUT fill factor = 2.0

-------------------------
Technical Discussion
-------------------------

By increasing this factor, the preconditioner becomes more accurate because more
terms in the preconditioner (pseudo-inverse) are retained. A value of 1.0 indicates that
the preconditioner would contain approximately the same number of nonzero entries as
the original matrix.

The two main parameters when using the ILUT preconditioner are this card and the
*Matrix drop tolerance* card. If the *Matrix drop tolerance* is 0.0, then this card
determines the size of the preconditioner. If *Matrix drop tolerance* is greater than 0.0,
then the approximate factorization is first created subject to this card’s restriction, and
then the drop tolerance is applied. This can result in a preconditioner with significantly
fewer nonzero entries.



