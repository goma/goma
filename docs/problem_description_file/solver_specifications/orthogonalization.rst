*********************
Orthogonalization
*********************

::

	Orthogonalization = {classic | modified}

-----------------------
Description / Usage
-----------------------

This optional card selects the orthogonalization scheme used internally for the gmres
solution algorithm (see the *Solution Algorithm* card). Valid options are

classic | classical
    Two steps of classical Gram-Schmidt orthogonalization.
modified
    A modified Gram-Schmidt orthogonalization.

If the *Orthogonalization* card is omitted, then the default selection is **classic.** Goma’s
parser will accept **classical** as equivalent to **classic.**

------------
Examples
------------

Following is a sample card:
::

	Orthogonalization = modified

-------------------------
Technical Discussion
-------------------------

By specifying **modified,** the user is greatly speeding up the **gmres** algorithm at the
expense of possibly losing convergence. A good indication that you should not have
used the **modified** setting is a premature “leveling off” of the sequence of residuals
produced internally within **gmres.**
