**************************
**Energy Weight Function**
**************************

::

   Energy Weight Function = {GALERKIN | SUPG} <float>

-----------------------
**Description / Usage**
-----------------------

This card specifies the weight function to be used on the weighted residual of the
energy equations. For high Peclet number cases, you may want to use a Petrov-
Galerkin formulation rather than a Galerkin formulation. Definitions of the input
parameters are as follows:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**GALERKIN**     |Name of the model for the weight functions for a full Galerkin formulation. This is the default when this   |
|                 |card is absent.                                                                                             |
|                 |                                                                                                            |
|                 | * <float> - the value of the weight function, a number between 0. and 1.; a value of 0. corresponds to     |
|                 |   **GALERKIN**.                                                                                            |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**SUPG**         |Name of the model for the weight functions for a streamwise upwinded Petrov-Galerkin formulation.           |
|                 |                                                                                                            |
|                 | * <float> - the value of the weight function, a number between 0. and 1.; a value of 1. corresponds to a   |
|                 |   full **SUPG**.                                                                                           |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:

::

   Energy Weight Function = GALERKIN 0.0

-------------------------
**Technical Discussion**
-------------------------

The **SUPG** weighting is applied only to the advective term in the Energy conservation
equation and Jacobian assembly.



--------------
**References**
--------------

No References.