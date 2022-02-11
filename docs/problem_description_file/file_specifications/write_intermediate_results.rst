******************************
Write Intermediate Results
******************************

::

	Write Intermediate Results = {yes | no}

-----------------------
Description / Usage
-----------------------

This optional card controls the output of intermediate results. The permissible values
for this card are

yes
    The code will output the latest Newton iteration to a file named ‘tmp.i.d’,
    where i is the Newton iteration number. The format of tmp.i.d will be
    similar to the ASCII results data described for the *GUESS file* and *SOLN
    file* cards. Also, the output EXODUS II database (see the *Output EXODUS II
    file* card) will accumulate the intermediate iterations as time planes of
    the solution.

no
    No intermediate results are written; only the last Newton iteration is
    written to the file named in the *SOLN file* card, and only the final
    converged iteration is output to the EXODUS II file.

------------
Examples
------------

Following is a sample card:
::

	Write Intermediate Results = no

-------------------------
Technical Discussion
-------------------------

This file is useful to guard against machine crashes or accidental job kills, particularly
for very large problems, as it can be used to restart a simulation (by using this file as the
*Guess file*). The intermediate results in the output EXODUS II database can be a useful
debugging tool, giving the analyst the ability to use highly relaxed Newton iterations to
see how a free boundary problem diverges.

