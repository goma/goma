******************
Number of Matrices
******************

::

	Number of Matrices = <integer>

-----------------------
**Description / Usage**
-----------------------

This card specifies the Number of Matrices for the current material. Note when running in
parallel all matrices should be on all materials

<integer>
   The number of matrices on this material

------------
**Examples**
------------


::


   ...
     Number of bulk species = 0

   Number of Matrices = 2

   MATRIX = 1
     Number of EQ   = 3
     ...
   
   MATRIX = 2
     Number of EQ   = 2
     ...

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



