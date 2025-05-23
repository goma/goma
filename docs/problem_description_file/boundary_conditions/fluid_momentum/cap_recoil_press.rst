********************
**CAP_RECOIL_PRESS**
********************

::

	BC = CAP_RECOIL_PRESS SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MOMENTUM)**

This boundary condition calculates the surface recoil from an evaporating metal alloy
component or water.

There are seven values in the <float_list>; definitions of the input parameters are as
follows:

==================== ======================================================
**CAP_RECOIL_PRESS** Name of the boundary condition (<bc_name>).
**SS**               Type of boundary condition (<bc_type>), where **SS** 
                     denotes side set in the EXODUS II database.
<bc_id>              The boundary flag identifier, an integer associated with
                     <bc_type> that identifies the boundary location (side 
                     set in EXODUS II) in the problem domain.
<float1>             This float is currently disabled.
<float2>             This float is currently disabled.
<float3>             Temperature at which the metal alloy begins to boil.
<float4>             Liquidus temperature of metal alloy.
<float5>             Reference temperature.
<float6>             Conversion scale for pressure.
<float7>             Conversion scale for temperature.
==================== ======================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = CAP_RECOIL_PRESS SS 1   0.0 0.0 3000.0 1623.0 0.0 1.0 1.0

-------------------------
**Technical Discussion**
-------------------------

Currently this boundary condition has coefficients for only iron and water. Several
required pieces of information to use this boundary condition are not in final form, and the user can expect future changes and improvements. This boundary condition is
designed for use with *Q_LASER_WELD*.



--------------
**References**
--------------

No References.