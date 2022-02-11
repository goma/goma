******************
**GRAD_LUB_PRESS**
******************

::

	BC = GRAD_LUB_PRESS SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/R_LUBP)**

This boundary condition card applies *free* boundary condition, akin to Papanastasiou et al. (1992) for the fluid momentum, at the boundary of a shell-element sheet. The boundary condition is applied to a sideset.

============== =======================================================
GRAD_LUB_PRESS Name of boundary condition.
SS             Type of boundary condition (<bc_type>), where SS
               denotes side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (node
               set in EXODUS II) in the problem domain.
<float1>       Flowrate in L^2/t. Usually set for NOBC effect.
============== =======================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = GRAD_LUB_PRESS SS   100 0.

This condition applied at sideset 100.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



