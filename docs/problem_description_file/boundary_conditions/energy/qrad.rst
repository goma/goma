********
**QRAD**
********

::

	BC = QRAD SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card specifies heat flux using both convective and radiative
terms. The <float_list> has four parameters; definitions of the input parameters are as
follows:

============ ==================================================================
**QRAD**     Name of the boundary condition (<bc_name>).
**SS**       Type of boundary condition (<bc_type>), where **SS** denotes
             side set in the EXODUS II database.
<bc_id>      The boundary flag identifier, an integer associated with
             <bc_type> that identifies the boundary location (side set in
             EXODUS II) in the problem domain.
<float1>     *h*, convective heat transfer coefficient.
<float2>     :math:`T_s`, sink temperature.
<float3>     :math:`\varepsilon`, total hemispherical emissivity.
<float4>     :math:`\sigma`, Stefan-Boltzmann constant.
============ ==================================================================

------------
**Examples**
------------

The following is a sample card:
::

   BC = QRAD SS 100 10.0 273.0 0.3 5.6697e-8

-------------------------
**Technical Discussion**
-------------------------

The heat flux definition for this card is a combined convective and radiative
formulation:

.. figure:: /figures/132_goma_physics.png
	:align: center
	:width: 90%

where *h* and :math:`T_s` are the convective heat transfer coefficient and the sink temperature,
and :math:`\varepsilon` and :math:`\sigma` are the total hemispherical emissivity and Stefan-Boltzmann constant,
respectively. The latter constant has been made an input parameter rather than a code
constant so that the user can specify its value in units that are consistent for the problem
being modeled.

The *QRAD* boundary condition can be used in place of *QCONV* by simply setting the
emissivity value to zero.




.. TODO -Line 48 has a picture that needs to be changed out with the correct equation.