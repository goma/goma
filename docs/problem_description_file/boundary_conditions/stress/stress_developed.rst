****************
STRESS_DEVELOPED
****************

::

	BC = STRESS_DEVELOPED SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(SIC/POLYMER_STRESS)**

This is a fully developed flow condition from Xie and Pasquali for viscoelastic flows

<bc_id>                                                          
   The boundary flag identifier, an integer associated  
   with <bc_type> that identifies the boundary location 
   (side set in EXODUS II) in the problem domain.       

------------
**Examples**
------------

::

   BC = STRESS_DEVELOPED SS 1

::

   BC = S11_1 NS 7   4.0   1.0


-------------------------
**Technical Discussion**
-------------------------

Replaces boundary contributions of stress equations with stress equations where the advective term is removed (:math:`v\cdot\nabla S \rightarrow 0`)


--------------
**References**
--------------

Xie, Xueying, and Matteo Pasquali. "A new, convenient way of imposing open-flow boundary conditions in two-and three-dimensional viscoelastic flows." Journal of non-newtonian fluid mechanics 122, no. 1-3 (2004): 159-176.