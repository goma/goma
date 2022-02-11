*********
**QCONV**
*********

::

	BC = QCONV SS <bc_id> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card specifies convective heat flux. Definitions of the input
parameters are as follows:

========= ===================================================================
**QCONV** Name of the boundary condition (<bc_name>).
**SS**    Type of boundary condition (<bc_type>), where **SS** denotes
          side set in the EXODUS II database.
<bc_id>   The boundary flag identifier, an integer associated with
          <bc_type> that identifies the boundary location (side set in
          EXODUS II) in the problem domain.
<float1>  *h*, heat transfer coefficient.
<float2>  :math:`T_s`, sink temperature.
========= ===================================================================

------------
**Examples**
------------

The following is a sample card:
::

   BC = QCONV SS 100 10.0 293.0

-------------------------
**Technical Discussion**
-------------------------

The convective heat flux is defined as

.. figure:: /figures/131_goma_physics.png
	:align: center
	:width: 90%

where *h* and :math:`T_s` are the convective heat transfer coefficient and the sink temperature, respectively.




.. TODO -Line 44 has a picture that needs to be changed out with the correct equation.