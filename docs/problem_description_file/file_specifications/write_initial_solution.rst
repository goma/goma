**************************
Write Initial Solution
**************************

::

	Write initial solution = {yes | no}

-----------------------
Description / Usage
-----------------------

This optional card controls the output of an initial solution prior to the start of a time
dependent simulation. The permissible values for this card are:

yes
    This value sets the flag WRITE_INITIAL_SOLUTION variable to “TRUE”. The
    initial solution vector will be written to an EXODUS II file and to an
    ASCII file (if the number of processors is not greater than
    DP_PROC_PRINT_LIMIT, currently set to 4 in *rf_io.h*).

no
    No initial solution is written.

------------
Examples
------------

Following is a sample card:
::

	Write Initial Solution = yes

-------------------------
Technical Discussion
-------------------------

This option is useful to activate when help is desired in debugging the startup portion of
a transient simulation.

