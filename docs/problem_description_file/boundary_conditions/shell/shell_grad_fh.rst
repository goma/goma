*****************
**SHELL_GRAD_FH**
*****************

::

	BC = BC = SHELL_GRAD_FH SS <bc_id> <float_list>

-----------------------
**Description / Usage**
-----------------------

**(SIC/R_SHELL_GRAD_FH)**

This boundary condition card sets a slope to the liquid film at the boundary of a shellelement sheet. The corresponding equation is EQ=shell_filmh. The boundary condition is applied to a node set.

=============== =====================================================
SHELL_GRAD_FH   Name of boundary condition.
SS              Type of boundary condition (<bc_type>), where SS
                denotes side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (node
                set in EXODUS II) in the problem domain.
<float1>        slope
=============== =====================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SSHELL_GRAD_FH SS   100 0.0

This condition applies a film slope of 0.0 at nodeset 100.

-------------------------
**Technical Discussion**
-------------------------

* he actual weighted residual equation that is applied to node on the surface 
  is

.. figure:: /figures/250_goma_physics.png
	:align: center
	:width: 90%

where :math:`\phi_i` is the finite element trial function, **n** is the outward-pointing normal to
the surface, :math:`\Sigma` and is the slope specified in the <float1>.

* The condition replaces the residual equation shell_filmh at the boundary.




.. TODO - Line 45 has an image that needs to be replaced with the equation. 