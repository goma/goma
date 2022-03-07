***********************
**Number of Materials**
***********************

::

	Number of Materials = <integer>

-----------------------
**Description / Usage**
-----------------------

This required card denotes how many material sections are contained in the *Problem Description File*. Each material section will have its own problem description, consisting of the following: *MAT* card, *Coordinate System* card, *Mesh Motion* card, *Number of bulk species* card, *Number of EQ* card, and zero or more equation cards. The input parameter is defined as

=========== ======================================================
<integer>   The number of MAT cards (i.e., material sections) that
            follow; this number must be greater than zero.
=========== ======================================================

If there are more MAT cards than specified by <integer>, *Goma* ignores all extras (i.e., the first *Number of Materials* material sections are read). If <integer> is set to -1, *Goma* automatically counts the MAT cards between the *Number of Materials* card and the *END OF MAT* card.

------------
**Examples**
------------

Following is a sample card, indicating that there are two materials:
::

   Number of Materials = 2

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.