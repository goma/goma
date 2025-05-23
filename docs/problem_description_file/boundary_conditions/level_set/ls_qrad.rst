***********
**LS_QRAD**
***********

::

	BC = LS_QRAD LS <integer> <float1> <float2> <float3> <float4>

-----------------------
**Description / Usage**
-----------------------

**(EMB/ENERGY)**

This boundary condition card specifies heat flux using both convective and radiative
terms. This heat flux value is applied as an “embedded” source term on the heat
conservation equation at the zero level set contour (cf. *BC = QRAD* for ALE surfaces).
It can be used both when subgrid or subelement integration is being used. The
<float_list> has four parameters; definitions of the input parameters are as follows:

A description of the input parameters follows:

============ ==============================================================
**LS_QRAD**  Name of the boundary condition.
**LS**       This string is used to indicated that this is a “boundary”
             condition is applied at an internal phase boundary defined
             by the zero contour of the level set function.
<integer>    An integer parameter than is permitted to take one of three
             values -1, 0, or 1. Depending upon the choice of this
             parameter the heat flux value is applied to the negative
             phase, both phase, or the positive phase, respectively.
             Details are given below.
<float1>     *h*, convective heat transfer coefficient.
<float2>     :math:`T_s`, sink temperature.
<float3>     :math:`\varepsilon`, total hemispherical emissivity.
<float4>     :math:`\sigma`, Stefan-Boltzmann constant.
============ ==============================================================

------------
**Examples**
------------

An example:
::

   BC = LS_QRAD LS 0 10.0 273.0 0.3 5.6697e-8

-------------------------
**Technical Discussion**
-------------------------

This is the level-set counterpart to *BC = QRAD* which is the same boundary condition
applied to a parameterized mesh surface. Please see the discussion of that input record
for the functional form of this boundary condition.



--------------
**References**
--------------

No References. 
