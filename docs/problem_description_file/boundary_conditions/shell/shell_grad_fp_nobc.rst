**********************
**SHELL_GRAD_FP_NOBC**
**********************

::

	BC = BC = SHELL_GRAD_FP_NOBC SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(WIC/R_SHELL_GRAD_FP_NOBC)**

This boundary condition card applies *free* boundary condition, akin to Papanastasiou et al. (1992) for the fluid momentum, at the boundary of a shell-element sheet. The boundary condition is applied to a sideset.

================== ======================================================
SHELL_GRAD_FP_NOBC Name of boundary condition.
SS                 Type of boundary condition (<bc_type>), where SS
                   denotes side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (node
                   set in EXODUS II) in the problem domain.
================== ======================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SSHELL_GRAD_FP_NOBC SS   100

This condition applied at sideset 100.

-------------------------
**Technical Discussion**
-------------------------

* The finite element formulation of the first equation of the film profile  
  equation boundary integral in the form of

.. figure:: /figures/249_goma_physics.png
	:align: center
	:width: 90%

* This condition is similar to the SHELL_GRAD_FP boundary condition, except
  that the condition is now a weak integrated condition that is *added* to the residual equations, instead of replacing them and the flux is no longer specified.



