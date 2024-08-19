*********************
Use AutoDiff Assembly
*********************

::

	Use AutoDiff Assembly = {yes | no}

-----------------------
Description / Usage
-----------------------

yes
    AutoDiff is used to compute the Jacobian's for supported BCs and equations. Currently this is quite limited but can be combined with non-autodiff equations.
no
    AutoDiff is not used

Default: no

------------
Examples
------------

Following is a sample card:
::

	Use AutoDiff Assembly = yes

-------------------------
Technical Discussion
-------------------------

Some formulations are quite complicated to compute Jacobian's by hand so only have AutoDiff routines. Mostly this is for fully-stabilized viscoelastic flow using the `SQRT_CONF` formulation.

Requires Goma to be built with Sacado from Trilinos.

--------------
References
--------------
