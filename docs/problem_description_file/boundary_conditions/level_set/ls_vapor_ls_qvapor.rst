**********************
**LS_VAPOR/LS_QVAPOR**
**********************

::

	BC = LS_VAPOR LS <integer> <float1> <float2> <float3> <float4>

-----------------------
**Description / Usage**
-----------------------

**(EMB/ENERGY)**

This boundary condition card specifies heat flux model derived from a laser welding
application. This particular contribution accounts for the energy lost by vapor flux.
This heat flux value is applied as an “embedded” source term on the heat conservation
equation at the zero level set contour (cf. *BC = Q_LASER_WELD* for ALE surfaces).
It can be used both when subgrid or subelement integration is being used. The
<float_list> has four parameters; definitions of the input parameters are as follows:

A description of the input parameters follows:

============= =========================================================
**LS_VAPOR**  Name of the boundary condition.
**LS**        This string is used to indicated that this is a “boundary”
              condition is applied at an internal phase boundary defined
              by the zero contour of the level set function.
<integer>     An integer parameter than is permitted to take one of three
              values -1, 0, or 1. Depending upon the choice of this
              parameter the heat flux value is applied to the negative
              phase, both phase, or the positive phase, respectively.
              Details are given below.
<float1>      T_scale. Temperature scaling.
<float2>      q_scale. Heat flux scaling.
============= =========================================================

------------
**Examples**
------------

An example:
::

   BC = LS_VAPOR LS 0 273. 1.

-------------------------
**Technical Discussion**
-------------------------

Currently this BC is hardwired to parameters (viz. heat capacitance, etc.) for iron. The
melting point temperature is taken from the material property “Liquidus Temperature”.
This boundary condition is still in the developmental stage. In using it is advisable to
be working with the Sandia Goma code team.



--------------
**References**
--------------

No References. 
