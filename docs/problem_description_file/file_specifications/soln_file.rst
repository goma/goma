*************
SOLN File
*************

::

	SOLN file = <file_name>

-----------------------
Description / Usage
-----------------------

This required card identifies the ASCII output file that will provide the initial guess for
continuation or time integration, where

<file_name>    
    Specifies the name of the output file, or if no file is desired, a value of
    **no** or **none** should be entered.

The current format of this ASCII file is a list of unformatted floating point numbers that
includes every degree of freedom in the problem in the order specified in the unknown
map. Other information (residual for that degree of freedom) may appear beyond the
first column of numbers in this file that is sometimes useful in determining the name
and location of the corresponding degree-of-freedom. If **no** or **none** is used in place of
the file name, no ASCII information is written.

------------
Examples
------------

Following is a sample card:
::

	SOLN file = soln.dat

-------------------------
Technical Discussion
-------------------------

This file represents the primary ASCII output of the *Goma* solution vector and the
primary way to continue or restart a solution from an ASCII file. (See *Write
Intermediate Solutions* for related information.) When a continuation run is performed,
this file is copied into the file specified in the *GUESS file* input card.

