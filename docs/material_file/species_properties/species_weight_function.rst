***************************
**Species Weight Function**
***************************

::

   Species Weight Function = {model_name} <float>

-----------------------
**Description / Usage**
-----------------------

This optional card is used to specify the weight functions to be used on the weighted
residual of the species convective diffusion equations. For high Peclet number cases,
you may want to use a Petrov-Galerkin formulation rather than a Galerkin formulation.

+-------------------+-------------------------------------------------------------------------------------+
|{model_name}       |Name of the formulation model. Valid entries are **GALERKIN**, for a full Galerkin   |
|                   |formulation, **SUPG**, for a streamwise upwinded Petrov-Galerkin formulation.        |
|                   |                                                                                     |
|                   | * <float> - the weight function parameter, chosen between 0. and 1.. The value 0.   |
|                   |   corresponds to **GALERKIN** weighting and 1. corresponds to a full **SUPG**.      |
+-------------------+-------------------------------------------------------------------------------------+

When this card is absent, the default {model_name} is **GALERKIN**.

------------
**Examples**
------------

The following is a sample input card:

::

   Species Weight Function = SUPG 0.5

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.