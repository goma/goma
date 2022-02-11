***********
**CURRENT**
***********

::

	BC = CURRENT SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

This card specifies the electrical current density at a given boundary.

Definitions of the input parameters are as follows:

============= ===============================================================
**CURRENT**   Name of the boundary condition (<bc_name>).
**SS**        Type of boundary condition (<bc_type>), where **SS** denotes
              side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (side set in
              EXODUS II) in the problem domain.
<float>       Value of current density 
              (in A/ :math:`m^2` or A/ :math:`cm^2`, depending on
              units of length scale used in the problem).
============= ===============================================================

------------
**Examples**
------------

An example input card:
::

   BC = CURRENT SS 1   -0.05

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



