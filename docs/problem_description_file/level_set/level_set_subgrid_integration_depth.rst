***************************************
Level Set Subgrid Integration Depth
***************************************

::

	Level Set Subgrid Integration Depth = <integer1>

-----------------------
Description / Usage
-----------------------

Subgrid integration is used to improve integration accuracy for all functions which
invoke a diffuse level-set interface representation of properties and surfaces. With
integration depths greater than zero the elements through which the zero level set
crosses are subdivided in a geometric way to achieve more accurate integration. Level-
1 depths implies the smallest grid size is 1/4 of the original, and a level-2 is 1/8th, and
so on. Please see usage nodes below.

<integer1>
    Level of integration depth. Default is zero. See usage notes.

------------
Examples
------------

This example sets the subgrid integration depth to two:
::

	Level Set Subgrid Integration Depth = 2

-------------------------
Technical Discussion
-------------------------

Each level of subgrid integration leads to precipitous growth in computational
load, especially in 3D. Level-2 seems to optimize accuracy and efficiency.
Levels higher than 2 is not recommended.

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
