**************
**VNORM_LEAK** 
**************

::

	BC = VNORM_LEAK SS <bc_id> <float1> <float2>

-----------------------
**Description / Usage**
-----------------------

This boundary condition card is used to specify a normal velocity boundary condition
with mass transfer on momentum equations. The flux quantity is specified on a per
mass basis so the heat and mass transfer coefficients are in units of L/t.

.. figure:: /figures/100_goma_physics.png
	:align: center
	:width: 90%

Definitions of the input parameters are as follows:

============== ================================================================
**VNORM_LEAK** Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS**
               denotes side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set
               in EXODUS II) in the problem domain.
<float1>       :math:`h_i`, mass transfer coefficient for bulk fluid 
               (n+ :math:`1^{th}` species).
<float2>       :math:`y^0_i`, driving force concentration in external phase.
============== ================================================================

------------
**Examples**
------------

The following is a sample input card:
::

     BC = VNORM_LEAK SS 1   1. 0.

-------------------------
**Technical Discussion**
-------------------------

This card is the equivalent of *KIN_LEAK* except it is solved for the normal component
of the momentum equation. Similar to *KIN_LEAK*, this flux provides an overall mass
transfer balance at an interface. Please refer to the technical discussion of *KIN_LEAK*
boundary card.




.. TODO - Line 17 has a photo that needs to be replaced with the real equation.