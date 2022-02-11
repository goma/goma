*************************
Matrix Drop Tolerance
*************************

::

	Matrix drop tolerance = <float>

-----------------------
Description / Usage
-----------------------

This optional card indicates to Aztec a drop tolerance to be used in conjunction with
preconditioners based on LU or on ILUT. The <float> input parameter is:

<float>
    **tol,** a floating point number ( :math:`tol \geq 0` ) that specifies the drop
    tolerance.

If the *Matrix drop tolerance* card is omitted, the default is **0.0.**

------------
Examples
------------

Following is a sample card:
::

	Matrix drop tolerance = 0.01

-------------------------
Technical Discussion
-------------------------

When constructing the partial factorization(s), any value less than tol is dropped. If set
to 0.0, then other parameters will govern preconditioner size and components (e.g.,
*Matrix ILUT fill factor* for the ILUT preconditioner).

The two main parameters when using the ILUT preconditioner are this card and the
*Matrix ILUT fill factor* card. The restrictions in *Matrix ILUT fill factor* take precedence
over the dropped entries caused by this card.



