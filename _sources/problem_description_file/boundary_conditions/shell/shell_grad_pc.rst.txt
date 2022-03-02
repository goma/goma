*****************
**SHELL_GRAD_PC**
*****************

::

	BC = BC = SHELL_GRAD_PC SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(WIC/R_SHELL_GRAD_PC)**

This boundary condition card allows the user to set volumetric flux of particles inside liquid film at the boundary of a shell-element sheet. The corresponding equation is EQ=shell_partc. The boundary condition is applied to a side set.

============= =========================================================
SHELL_GRAD_PC Name of boundary condition.
SS            Type of boundary condition (<bc_type>), where SS
              denotes side set in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (node
              set in EXODUS II) in the problem domain.
<float1>      particle flux
============= =========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SHELL_GRAD_PC SS   100 1.

This condition applied at sideset 100. and sets a particle flux to 1.0

-------------------------
**Technical Discussion**
-------------------------

* The actual weighted residual equation that is applied to node on the surface 
  is

.. figure:: /figures/252_goma_physics.png
	:align: center
	:width: 90%

where :math:`\phi_i` is the finite element trial function, **n** is the outward-pointing normal to the surface, :math:`J_p` and is the particles flux specified in the <float1>.

* The condition replaces the residual equation shell_partc at the boundary.




.. TODO - Line 45 has an image that needs to be replaced with the equation. 
