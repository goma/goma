**********************
**SHELL_GRAD_FH_NOBC**
**********************

::

	BC = BC = SHELL_GRAD_FH_NOBC SS <bc_id>

-----------------------
**Description / Usage**
-----------------------

**(WIC/R_SHELL_GRAD_FH_NOBC)**

This boundary condition card applies *free* boundary condition, akin to Papanastasiou et al. (1992) for the fluid momentum, at the boundary of a shell-element sheet, in terms of the slope of a thin Reynolds film. The boundary condition is applied to a sideset.

================== =========================================================
SHELL_GRAD_FH_NOBC Name of boundary condition.
SS                 Type of boundary condition (<bc_type>), where SS
                   denotes side set in the EXODUS II database.
<bc_id>            The boundary flag identifier, an integer associated with
                   <bc_type> that identifies the boundary location (node
                   set in EXODUS II) in the problem domain.
================== =========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SSHELL_GRAD_FH_NOBC SS   100

This condition applied at sideset 100.

-------------------------
**Technical Discussion**
-------------------------

* The finite element formulation of the second equation of film profile 
  equation generates boundary integral in the form of

.. figure:: /figures/251_goma_physics.png
	:align: center
	:width: 90%

* This condition is similar to the SHELL_GRAD_FH boundary condition, except
  that the condition is now a weak integrated condition that is *added* to the 
  residual equations, instead of replacing them and the flux is no longer 
  specified.




.. TODO - Line 44 has an image that needs to be replaced with the equation. 
