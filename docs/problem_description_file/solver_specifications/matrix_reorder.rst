******************
Matrix Reorder
******************

::

	Matrix reorder = {none | rcm}

-----------------------
Description / Usage
-----------------------

This optional card determines whether RCM (Reverse Cuthill-McKee) reordering of
the linear system is to be performed. Valid options are:

none
    the equations are not reordered.
rcm
    the equations are reordered using an RCM scheme.

If the *Matrix reorder* card is omitted, then the default selection is **none.**

------------
Examples
------------

Following is a sample card:
::

	Matrix reorder = rcm

-------------------------
Technical Discussion
-------------------------

Note that reordering frequently is helpful in achieving convergence for iterative
solution of linear systems. In a few instances, however, *Goma* users have noted that
RCM reordering hinders convergence for selected problems. The default for *Goma* is to
not use the RCM reordering so that quantitatively comparable results are obtained
using either Aztec 1 (which did not have RCM reordering as an option) or Aztec 2.x. In
summary, users are encouraged to try RCM reordering when using iterative solvers,
foregoing the option only as a further resort in the face of repeated convergence
failures.




