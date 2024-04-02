****************
Print 3D BC Dup
****************

::

	Print 3D BC Dup = <integer>

-----------------------
Description / Usage
-----------------------

This optional card specifies the level of information output to the 3D BC dup file "bc3D_output.txt" 

.. tabularcolumns:: |l|L|

==============  ===============================================================
Level           Results Output
==============  ===============================================================
0               No output (default).
1               BC Duplication information
2               Same as level 1 with additional rotation and node Duplication
                debugging
==============  ===============================================================


------------
Examples
------------

Following is a sample card:
::

	Print 3D BC Dup = 1


-------------------------
Technical Discussion
-------------------------