*****************************
**Number of Phase Functions**
*****************************

::

	Number of phase functions = {integer}

-----------------------
**Description / Usage**
-----------------------

Activates generalized phase function capability. Currently, the number of phase
functions cannot exceed five. Phase function fields are essentially identical to level set
fields, but more than one can be activated for various purposes. Please see technical
discussion below.

------------
**Examples**
------------

A sample input card is:
::

	Number of phase functions = 1

-------------------------
**Technical Discussion**
-------------------------

Various uses of the phase function approach have been explored. To track multiple
interface types from multiple fluids requires more than one level-set field. This
capability can also be deployed for tracking imprinted solid surfaces (moving) together
with capillary free surfaces. Consult the tutorials.


--------------
**References**
--------------

GT-026.3 GOMA’s Overset Mesh Method: User Tutorial, November 19 2003. P. R.
Schunk and E. D. Wilkes
