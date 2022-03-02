****************************
Level Set Adaptive Order
****************************

::

	Level Set Adaptive Integration = <integer1>

-----------------------
Description / Usage
-----------------------

To be used with Subelement adaptive integration to improve integration accuracy. Does
not work with subgrid integration or basic level-set. Requires a sharp interface, viz.
level-set length scale of zero. Please see usage nodes below.

<integer1>
    Adaptive integration order. Single positive integer greater than zero.
    Default value is 3.

------------
Examples
------------

This example invokes the subelement adaptive integration order:
::

	Level Set Adaptive Integration = YES

::

	Level Set Adaptive Order = 2

--------------
References
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
