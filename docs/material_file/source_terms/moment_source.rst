***************
Moment Source
***************

::

   Moment Source = {model_name} <float_list> [varies]

-----------------------
Description / Usage
-----------------------

This required card is used to specify the model for the source term on the energy
equation. Definitions of the input parameters are as follows:

{model_name}
    Name of the model for the source term on the energy equation. The permissible values
    are
    * FOAM_PMDI_10
    * CONSTANT_GROWTH
    * FOAM_PBE

<float_list> 
    One or more floating point numbers (<float1> through <floatn>) whose values are
    determined by the selection for {model_name}. Note that not all models have a
    <float_list>


Source-term model choices and their parameters are discussed below. WARNING:
make sure the equation term multipliers for the source terms are set to unity (see the
Equation Cards segment in the previous chapter).

FOAM_PMDI_10
    Foam PMDI_10 model from Ortiz et al.
    * <float1> Growth rate coefficient, multiplicative with Moment Growth Rate Kernel
    * <float2> Coalescence kernel (Unused, set to 1.)

CONSTANT_GROWTH
    Constant growth rate
    * <float1> Growth rate coefficient

------------
**Examples**
------------

The following is a sample input card:

::

   Moment Source = FOAM_PMDI_10 1. 1.

-------------------------
**Technical Discussion**
-------------------------


--------------
**References**
--------------

Ortiz, Weston, et al. "Population balance modeling of polyurethane foam formation with pressure‐dependent growth kernel." AIChE Journal (2021): e17529.