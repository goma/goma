***********************
Polymer Weight Function
***********************

::

   Polymer Weight Function = {GALERKIN | SUPG}

-----------------------
**Description / Usage**
-----------------------

This optional card is used to specify the weight function for the polymer stress
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

The following is a sample card that set the polymer weight function to **SUPG** and
demonstrates the required cards.

::

   Polymer Weight Function = SUPG

::

   Polymer Weighting = CONSTANT 0.1

The following is a sample card that set the polymer weight function to **GALERKIN**.

::

   Polymer Weight Function = GALERKIN

-------------------------
**Technical Discussion**
-------------------------

No Discussion.




