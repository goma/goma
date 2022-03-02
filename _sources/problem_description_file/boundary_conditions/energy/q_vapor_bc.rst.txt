**************
**Q_VAPOR_BC**
**************

::

	BC = Q_VAPOR SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used to specify heat loss due to evaporation. It is
typically used in conjunction with Q_LASER_WELD.

=========== ===============================================================
**Q_VAPOR** Name of the boundary condition (<bc_name>).
**SS**      Type of boundary condition (<bc_type>), where **SS**
            denotes side set in the EXODUS II database.
<bc_id>     The boundary flag identifier, an integer associated with
            <bc_type> that identifies the boundary location (side set
            in EXODUS II) in the problem domain.
<float1>    Temperature scale.
<float2>    Energy unit scale.
=========== ===============================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = Q_VAPOR SS 10 100. 10.

-------------------------
**Technical Discussion**
-------------------------

This condition is turned on above the boiling point, which is story in the melting 
point solidus temperature.




