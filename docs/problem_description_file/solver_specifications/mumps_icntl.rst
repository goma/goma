***********
MUMPS ICNTL
***********

::

	MUMPS ICNTL <int>  = <value>

-----------------------
Description / Usage
-----------------------

The MUMPS ICNTL parameter is used to set the control parameters for the MUMPS
solver. The parameters are set using the ICNTL array, which is a 1-based array
of integers. The first <int> of the card corresponds to the parameter number,
and the second <value> corresponds to the value to be set.


------------
Examples
------------

Following is a sample card:
::

    # Enable extra information from MUMPS solver
	MUMPS ICNTL 4 = 2

