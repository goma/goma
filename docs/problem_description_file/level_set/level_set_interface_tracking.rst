********************************
Level Set Interface Tracking
********************************

::

	Level Set Interface Tracking = {yes | no}

-----------------------
Description / Usage
-----------------------

Activates (or deactivates) embedded interface tracking by the level set method. When
activated, the set of cards specifying level set run parameters are read; these should
appear in the input deck following this card. Also when activated a “level_set”
equation type should be included in the list of equations identified in the equations
section.

------------
Examples
------------

A sample input card is:
::

	Level Set Interface Tracking = yes

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
