*****************
**SHELL_GRAD_FP**
*****************

::

	BC = BC = SHELL_GRAD_FP SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/R_SHELL_GRAD_FP)**

This boundary condition card applies a volumetric flux of liquid film to the boundary of a shell-element sheet. The corresponding equation is EQ=shell_filmp. The boundary condition is applied to a node set.

============= =========================================================
SHELL_GRAD_FP Name of boundary condition.
SS            Type of boundary condition (<bc_type>), where SS
              denotes side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (node
              set in EXODUS II) in the problem domain.
<float1>      volumetric flux
============= =========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SSHELL_GRAD_FP SS   100 0.0

This condition applies a particles volume flux of 0.0 at nodeset 100.

-------------------------
**Technical Discussion**
-------------------------

The actual weighted residual equation that is applied to node on the surface is

.. figure:: /figures/248_goma_physics.png
	:align: center
	:width: 90%

where :math:`\phi_i` is the finite element trial function, **n** is the outward-pointing normal to
the surface, and q is the volumetric flux specified in the <float1>. Careful attention should be given for the sign of q. The **sign convention** is that q is **positive** when the flow is **exiting** the boundary and **negative** when **entering** the boundary.

The condition replaces the residual equation shell_filmp at the boundary.



