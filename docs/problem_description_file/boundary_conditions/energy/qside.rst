*********
**QSIDE**
*********

::

	BC = QSIDE SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition card is used to specify a constant heat flux. Definitions of the
input parameters are as follows:

========== =================================================================
**QSIDE**  Name of the boundary condition (<bc_name>).
**SS**     Type of boundary condition (<bc_type>), where **SS** denotes
           side set in the EXODUS II database.
<bc_id>    The boundary flag identifier, an integer associated with
           <bc_type> that identifies the boundary location (side set in
           EXODUS II) in the problem domain.
<float1>   Value of heat flux. A positive value implies that energy is
           being added to system; a negative value implies energy is
           being taken from the system through the boundary.
========== =================================================================

------------
**Examples**
------------

The following is a sample card:
::

   BC = QSIDE SS 22   1.50

-------------------------
**Technical Discussion**
-------------------------

The mathematical form of the boundary condition.

.. figure:: /figures/133_goma_physics.png
	:align: center
	:width: 90%




.. TODO -Line 45 has a picture that needs to be changed out with the correct equation.