**************************
Moment Growth Kernel
**************************

::

   Moment Growth Kernel = {CONSTANT | VISCOSITY_SCALED | VISCOSITY_PRESSURE_SCALED} <float1> [float2]

-----------------------
Description / Usage
-----------------------

This card specifies the growth kernel to be used for the
Moment equations. Currently only used for `PMDI_10` QMOM models.
Definitions of the input parameters are as follows:

CONSTANT     
    Constant growth rate
    * <float1> the prefactor

VISCOSITY_SCALED
    Growth rate scaled by viscosity :math:`G = G_0 \eta_0 / \mu`, for use with `PMDI_10` models
    * <float1> - :math:`G_0`

VISCOSITY_PRESSURE_SCALED
    Growth rate scaled by viscosity and pressure :math:`G = G_0 (\eta_0 / \mu) / (p - p_ref)^2`, for use with `PMDI_10` models

    * <float1> :math:`G_0`
    * <float2> :math:`p_ref`
 
------------
Examples
------------

The following is a sample input card:

::

   Moment Growth Kernel = CONSTANT 1e-6

   Moment Growth Kernel = VISCOSITY_SCALE 1e-3

   Moment Growth Kernel = VISCOSITY_PRESSURE_SCALED 1e-8

-------------------------
Technical Discussion
-------------------------

Note this is multiplicative with the growth rate in the `FOAM_PMDI_10` moment source model

--------------
References
--------------

Ortiz, Weston, et al. "Population balance modeling of polyurethane foam formation with pressure‐dependent growth kernel." AIChE Journal (2021): e17529.