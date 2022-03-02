****************
**Output Level**
****************

::

	Output Level = <integer>

-----------------------
**Description / Usage**
-----------------------

This optional card specifies the level of diagnostic information output to the file
stderr. The permissible values for <integer> are **0** through **4**, depending on the level
of informational (debugging) output desired; higher values of the output level will
produce more diagnostic information on the stdout and stderr output channels.
The default output level is **0**. Specific output is summarized below.

.. tabularcolumns:: |l|L|

==============  ===============================================================
Level           Results Output
==============  ===============================================================
0               No diagnostic output (default).
1               Identifies the degree of freedom, the solution variable, and
                node at which the maximum value of norm is present.
2,3,4           Currently unused; available for developer output
                specification.
==============  ===============================================================

------------
**Examples**
------------

Following is a sample card:
::

	Output Level = 1

-------------------------
**Technical Discussion**
-------------------------

This specification allows the developer a means to output specific information that
would be helpful in diagnosing problems in the software. Currently, the output options
are limited.

