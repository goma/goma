****************************************
**SHELL_OPEN_PRESS, SHELL_OPEN_PRESS_2**
****************************************

::

	BC = SHELL_OPEN_PRESS NS <bc_id> <float_list>

::

	BC = SHELL_OPEN_PRESS_2 NS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(DC/R_SHELL_SAT_OPEN or SHELL_SAT_OPEN_2)**

This Dirichlet boundary condition card applies a shell liquid phase pressure to the boundary of a shell-element sheet. The corresponding equation is EQ=shell_sat_open or correspondingly shell_sat_open_2, depending on which layer. The boundary condition is applied to a node set.

================ ===================================================
SHELL_OPEN_PRESS Name of boundary condition.
NS               Type of boundary condition (<bc_type>), where NS
                 denotes node set in the EXODUS II database.
<bc_id>          The boundary flag identifier, an integer associated with
                 <bc_type> that identifies the boundary location (node
                 set in EXODUS II) in the problem domain.
<float1>         SHELL_OPEN_PRESSURE, the value of the liquid
                 phase pressure at the boundary.
[float2]         Optional floating point number set between 0.0 and 1.0
                 which serves as a flag to the code for a Dirichlet
                 boundary condition. If this value is present, and is not
                 1.0, the condition is applied as a residual equation.
                 Otherwise, it is “hard-set” condition and is eliminated
                 from the matrix. The residual method must beused
                 when this Dirichlet boundary condition is used as a
                 parameter in automatic continuation sequences.
================ ===================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SHELL_OPEN_PRESS 100 1.0

This boundary condition is applied at nodeset 100.

-------------------------
**Technical Discussion**
-------------------------

* The equation applied at the specified nodeset in place of the shell-sat-open
  equation.



--------------
**References**
--------------

No References.