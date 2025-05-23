*****************************
Matrix Residual Norm Type
*****************************

::

	Matrix residual norm type = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card selects the type of norm that is used to measure the size of the
residuals occurring during the solution of the linear matrix system :math:`r \left( z \right) = b - Az`, where
:math:`z` is an approximation to the solution :math:`x` of the linear matrix problem :math:`Ax = b`. The types
of norms used by the linear solver are controlled by values of {char_string}:

.. tabularcolumns:: |l|L|

==================  ===============================================================================================================
{char_string}       Norm type
==================  ===============================================================================================================
**r0**              :math:`\frac{\lVert r \rVert_2}{\lVert r^0 \rVert_2}`
**rhs**             :math:`\frac{\lVert r \rVert_2}{\lVert b \rVert_2}`
**Anorm**           :math:`\frac{\lVert r \rVert_2}{\lVert A \rVert_{\infty}}`
**sol**             :math:`\frac{\lVert r \rVert_{\infty}}{\lVert A \rVert_{\infty} \lVert x \rVert_1 + \lVert b \rVert_{\infty} }`
**noscaled**        :math:`\lVert r \rVert_2`
==================  ===============================================================================================================

The (0) superscript for the **r0** specification indicates the initial value of the residual.

If the *Matrix residual norm type* card is omitted, the default is **r0.**

------------
**Examples**
------------

Following is a sample card:
::

	Matrix residual norm type = r0

-------------------------
**Technical Discussion**
-------------------------

For direct factorization linear solution algorithms, the norm should become very small
in the single iteration that is performed. This card is more pertinent when an iterative
solution algorithm has been specified.

Note the distinction between the residual for the overall global Newton iteration and
use of the term residual to describe an aspect of the linear solver iteration. For the linear
matrix systems, a residual :math:`r` may be computed for any guess of the solution to :math:`Ax = b` as
:math:`r(z) = b - Az`. If :math:`z = x`, the actual solution, then the residual is zero; otherwise, it is 
some vector with a nonzero norm.
