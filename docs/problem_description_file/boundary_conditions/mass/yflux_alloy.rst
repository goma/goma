***************
**YFLUX_ALLOY**
***************

::

	BC = YFLUX_ALLOY SS <bc_id> <integer1> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary condition card calculates the surface integral for a mass flux transfer
model for the evaporation rate of molten metal.

The <float_list> requires six values; definitions of the input parameters are as follows:

=============== =================================================================
**YFLUX_ALLOY** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS**
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set
                in EXODUS II) in the problem domain.
<integer1>      Species number.
<float1>        Liquidus temperature of metal alloy, :math:`T_m`.
<float2>        Base Concentration, :math:`y^\infty`.
<float3>        Coefficient :math:`c_0`.
<float4>        Coefficient :math:`c_1`.
<float5>        Coefficient :math:`c_2`.
<float6>        Coefficient :math:`c_3`.
=============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_ALLOY SS 10 0 1623.0 0.5 0.01 -1e-3 1e-4 -1e-5

-------------------------
**Technical Discussion**
-------------------------

Basically the difference between this model and the simple convective mass transfer
coefficient (say :math:`k_i` for *YFLUX*) is that the transfer coefficient here (the exponential term) has a cubic dependence on temperature.

.. figure:: /figures/153_goma_physics.png
	:align: center
	:width: 90%




.. TODO - Line 52 has a photo that needs to be replaces with the proper equation.