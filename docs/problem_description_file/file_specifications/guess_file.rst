**************
Guess File
**************

::

	GUESS file = <file_name>

-----------------------
Description / Usage
-----------------------

This required card identifies the input file that provides the initial guess for the solution
vector for continuation or time integration, where

<file_name>    
    Specifies the exact name of the file and can be any file name.

The file <file_name> is read by *Goma* only if the value of the *Initial Guess* (next
section on *General Specifications*) card is set to read. The current format of this ASCII
file is a list of unformatted floating point numbers (the solution variable followed by
the residual value for that degree of freedom) in the order of the unknown map; this is
the same format as the file described in the *SOLN* file card. A solution file from a
previous simulation may be used.

------------
Examples
------------

Following is a sample card:
::

	GUESS file = contin.dat

-------------------------
Technical Discussion
-------------------------

This file is typically a copy of the *SOLN file* thus being an exact replica of it. It
represents the only way to continue a previous solution from an ASCII file. Typically a
continuation proceeds from a converged solution but the result from an intermediate
solution could also be used; the user is cautioned about the potential difficulties of
restarting from non-converged solution. (See *Initial Guess* card about (re-)starting from
a binary file.)

