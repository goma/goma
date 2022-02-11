**********************
**LS_RECOIL_PRESSURE**
**********************

::

	BC = LS_RECOIL_PRESSURE LS <integer> <float1> <float2> <float3> <float4>

-----------------------
**Description / Usage**
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition card specifies heat flux model derived from a laser welding
application.. This heat flux value is applied as an “embedded” source term on the heat
conservation equation at the zero level set contour (cf. *BC = CAP_RECOIL_PRESS* for
ALE surfaces). It can be used both when subgrid or subelement integration is being
used. The <float_list> has seven parameters; definitions of the input parameters are as
follows:

A description of the input parameters follows:

====================== ============================================================
**LS_RECOIL_PRESSURE** Name of the boundary condition.
**LS**                 This string is used to indicated that this is a “boundary”
                       condition is applied at an internal phase boundary defined
                       by the zero contour of the level set function.
<integer>              An integer parameter than is permitted to take one of three
                       values -1, 0, or 1. Depending upon the choice of this
                       parameter the heat flux value is applied to the negative
                       phase, both phase, or the positive phase, respectively.
                       Details are given below.
<float1>               This float is currently disabled.
<float2>               This float is currently disabled.
<float3>               This float is currently disabled.
<float4>               Disabled. The boiling temperature is set to the melting
                       point of the solidus. Use the material property “Solidus
                       Temperature” card for this.
<float5>               This float is currently disabled.
<float6>               Conversion scale for pressure.
<float7>               Conversion scale for temperature.
====================== ============================================================

------------
**Examples**
------------

-------------------------
**Technical Discussion**
-------------------------

Currently this boundary condition has coefficients for only iron and water. Several
required pieces of information to use this boundary condition are not in final form, and
the user can expect future changes and improvements. This boundary condition is
designed for use with *LS_QLASER*.



--------------
**References**
--------------

No References. 
