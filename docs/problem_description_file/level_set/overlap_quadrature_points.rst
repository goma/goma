*****************************
Overlap Quadrature Points
*****************************

::

	Overlap Quadrature Points = <integer1>

-----------------------
Description / Usage
-----------------------

To be used with the overset grid capability. This function sets the number of overlap
quadrature points with this capability. See GT-026 for more details.

<integer1>
    Overlap quadrature points. Single positive integer greater than zero.
    Default value is 3.

------------
Examples
------------

This example invokes the number of overlapping quadrature points:
::

	Overlap Quadrature Points = 2

-------------------------
Technical Discussion
-------------------------

Please consult the overset grid capability tutorial for futher discussion. (Ref.
below). This is to be use with AC_OVERLAP, or the augmenting condition of
type AC = OV.

--------------
References
--------------

GT-026.4: “GOMA’s Overset Mesh Method”, P R. Schunk and E. D. Wilkes, 11 Jan.
2006
