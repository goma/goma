********************************
Spalart Allmaras Weight Function
********************************

::

   Spalart Allmaras Weight Function = {GALERKIN | SUPG}

-----------------------
**Description / Usage**
-----------------------

This optional card is used to specify the weight function for the Spalart Allmaras turbulent model
equation. Valid options are

+-----------------+----------------------------------------------------------------------------------------------------------------+
|**GALERKIN**     |Uses a Galerkin weight-function for the stress equation. This option is the default if this card is not present.|
+-----------------+----------------------------------------------------------------------------------------------------------------+
|**SUPG**         |Uses a streamline upwind Petrov-Galerkin weight-function for the stress equation. If this option is chosen, a   |
|                 |weight must be specified via the Polymer Weighting card.                                                        |
+-----------------+----------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------


::

   Spalart Allmaras Weight Function = SUPG 1.0

::

   Spalart Allmaras Weight Function = GALERKIN


-------------------------
**Technical Discussion**
-------------------------

No Discussion.




