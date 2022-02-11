********************************
Matrix Factorization Overlap
********************************

::

	Matrix factorization overlap = {char_string}

-----------------------
Description / Usage
-----------------------

This optional card determines how much matrix factorization overlap occurs with other
processors. This specification is only relevant for parallel computations. The valid
options for {char_string} are:

none
    No augmentation is performed, equivalent to a setting of k=0. This is the
    default.
diag
    Augment the processor’s local matrix to include the diagonal (MSR) or
    diagonal blocks (VBR) for external rows.
k
    Augment the processor’s local matrix to include external rows. The rows are
    selected by examining non-zero columns from the current local system that
    refer to offprocessor unknowns, and including the rows associated with
    those off-processor unknowns. This process is repeated *k* times, where
    k > 0. When complete, all non-zero columns whose associated rows have not
    been included are discarded. A value of 0 is equivalent to a setting of
    none.

If the *Matrix factorization overlap* card is omitted, the default is **none.**

------------
Examples
------------

Following is a sample card:
::

	Matrix factorization overlap = 1

-------------------------
Technical Discussion
-------------------------

This optional card determines how much a processor’s local matrix is to be augmented
with information from adjacent processors during the approximate factorizations used
to build preconditioners. This card should be omitted or given a value of **none** for serial
executions.



