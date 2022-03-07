**************************
Moment Coalescence Kernel
**************************

::

   Moment Coalescence Kernel = {ADDITION | VISCOSITY_SCALED_ADDITION | VISCOSITY_ADDITION_BUBBLE_RATIO} <float1>

-----------------------
Description / Usage
-----------------------

This card specifies the coalescence kernel to be used for the
Moment equations. Currently only used for `PMDI_10` QMOM models.
Definitions of the input parameters are as follows:

ADDITION     
    :math:`\beta(v,v') = \beta_0 (v + v')`
    * <float1> :math:`\beta_0`

VISCOSITY_SCALED_ADDITION
    :math:`\beta(v,v') = \beta_0 (\eta_0 / \mu) (v + v')`
    * <float1> :math:`\beta_0`

VISCOSITY_ADDITION_BUBBLE_RATIO
    :math:`\beta(v,v') = \beta_0 (\eta_0 / \mu) (v + v') (v / (v' + \epsilon))`
    * <float1> :math:`\beta_0`

 
------------
Examples
------------

The following is a sample input card:

::

   Moment Coalescence Kernel = ADDITION 1e-12

   Moment Coalescence Kernel = VISCOSITY_SCALED_ADDITION 1e-3

   Moment Coalescence Kernel = VISCOSITY_ADDITION_BUBBLE_RATIO 1e-8

-------------------------
Technical Discussion
-------------------------

Currently only for `FOAM_PMDI_10` moment source model

--------------
References
--------------

Ortiz, Weston, et al. "Population balance modeling of polyurethane foam formation with pressure‐dependent growth kernel." AIChE Journal (2021): e17529.