**********
MUMPS CNTL
**********

::

	MUMPS ICNTL <int>  = <value>

-----------------------
Description / Usage
-----------------------

The MUMPS CNTL parameter is used to set the control parameters for the MUMPS
solver. The parameters are set using the CNTL array, which is a 1-based array of
reals (double's in Goma's case). The first <int> of the card corresponds to the
parameter number, and the second <value> corresponds to the value to be set.
