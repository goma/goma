**************************
Decomposition Type
**************************

::

	Decomposition Type = {rcb | kway}

-----------------------
Description / Usage
-----------------------

This optional card controls which builtin METIS method for decomposition in parallel

rcb
    Recurisve Bisection, default when less than 8 processors and this card is not specified

kway
    KWAY default when 8 or more processors and this card is not specified

------------
Examples
------------

Following is a sample card:
::

	Decomposition Type = kway

-------------------------
Technical Discussion
-------------------------

Also available from the command line with :code:`-rcb, -kway`. 


