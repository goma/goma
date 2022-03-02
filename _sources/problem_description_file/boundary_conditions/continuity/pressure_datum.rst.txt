******************
**PRESSURE DATUM**
******************

::

	PRESSURE DATUM = <integer> <float>

-----------------------
**Description / Usage**
-----------------------

This card is used to set a hydrodynamic pressure datum on fluid mechanics problems
that contain no implicit or explicit boundary conditions on stress or pressure.
Definitions of the input parameters are as follows:

========== ========================================================
<integer>  Element number on which the datum is set. This number
           should correspond to that shown when viewing the mesh,
           less one, as the numbering convention in the C language
           starts at zero rather than at one.
<float>    Value of the hydrodynamic pressure datum.
========== ========================================================

Noteworthy is that this card is optional, and if used, is placed outside the BC section
and just below it.

------------
**Examples**
------------

Following is a sample card:
::

   PRESSURE DATUM = 10 1.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



