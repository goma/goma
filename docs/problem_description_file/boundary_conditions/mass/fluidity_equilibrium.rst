********************
FLUIDITY_EQUILIBRIUM
********************

::

	BC = FLUIDITY_EQUILIBRIUM SS <bc_id> <species index>

-------------------
Description / Usage
-------------------

**(SIC/MASS)**

Definitions of the input parameters are as follows:


FLUIDITY_EQUILIBRIUM
   Name of the boundary condition

SS
   This is applied to a side set

<bc_id>
   The SS id where we apply this BC

<species index>
   The 0-based species index associated with the FLUIDITY source term

------------
**Examples**
------------

Following are two sample cards:
::

   BC = FLUIDITY_EQUILIBRIUM SS   2

--------------------
Technical Discussion
--------------------

This replaces the boundary equations with the equivalent equation with
:math:`\mathbf{u}\cdot\nabla c = 0` replacing the advection term.



----------
References
----------
