*******************
**SH_FLUID_STRESS**
*******************

::

	BC = SH_FLUID_STRESS SS <bc_id> <float>

-----------------------
**Description / Usage**
-----------------------

**(PCC/VECTOR MOMENTUM)**

Used for fluid-structure interaction problems with structural shell elements, the *SH_FLUID_STRESS* condition equates the normal traction (the tangential and normal force components, per unit area) between adjacent fluid and solid structure. This condition is only to be used on boundaries between regions of ARBITRARY mesh motion with fluid momentum equations: see *Mesh Motion* and *EQ* cards. With this boundary condition, the local residual and Jacobian contributions from the fluid mechanics momentum equations (on the *ARBITRARY* side of the boundary) are added into weak form of the residual and Jacobian entries for the solid structural equations (*see EQ = shell_curvature and EQ = shell_tension*). All elements on both sides of the interface must have the same element type, i.e., the same order of interpolation and basis functions, e.g., Q1 fluid and Q1 (bar element) for shell. Q2 fluid momentum and Q2 (bar element) for the shell equations. Also, such interfaces must be defined as a
mesh side set attached to the bulk fluid elements (most mesh generators will not allow for side sets in bar or sheet elements).

Definitions of the input parameters are as follows:

=================== ========================================================
**SH_FLUID_STRESS** Name of the boundary condition (<bc_name>).
**SS**              Type of boundary condition (<bc_type>), where **SS** 
                    denotes side set in the EXODUS II database.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set in EXODUS II) in the problem domain.
<float>             Scale factor for stress balance for non-dimensionalization.
                    This parameter, which multiplies the liquid phase
                    contribution and should be set to 1.0 if there is no
                    nondimensional treatment.
=================== ========================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = SH_FLUID_STRESS SS 5 1.0

In this example, side set 5 is a boundary between a solid blade and a liquid; material 2 is the rubber blade, and material 1 is the fluid. Along that blade, a companion boundary condition of the form

::

   BC = NO_SLIP SS 5 2 1

should also be applied.

-------------------------
**Technical Discussion**
-------------------------

The functional form of the boundary condition is:

.. figure:: /figures/247_goma_physics.png
	:align: center
	:width: 90%

where :math:`\underline{T}` is the fluid phase stress tensor given by any one of the specified fluid-phase
constitutive equations, and :math:`\underline\sigma` is the solid-phase stress tensor, also given by any one of
the solid-phase constitutive equation (see material file specifications).
:math:`\lambda` is a scaling
factor that defaults to unity (and is usually best taken as such unless some scaling is
invoked).

This balance is applied to the weak form of the solid-phase momentum residuals, from
the fluid phase, viz. in the fluid-phase, the fluid-stress at the interface is added to the
solid-phase momentum residuals. As mentioned above, this condition usually needs to
be supplemented by a statement of mass conservation across the interface, which will
depend on whether the solid phase is of *CONTINUOUS* or *POROUS* media (see *Media Type* card).


--------
**FAQs**
--------

* **Troubleshooting 1:** This boundary condition requires that the side set 
  contain elements from both the fluid and the solid side of the interface. For the FASTQ tool, this is the default case; for CUBIT and possibly other related tools, this can be forced on the side set definition options. Interestingly, the boundary condition does work if the side set is attached to the fluid phase only, but just due to the way in which it is applied.

* **Troubleshooting 2:** This boundary condition does not enforce mass
  conservation. A combination of NO_SLIP or VELO_NORMAL/VELO_TANGENT must be invoked to achieve a material surface. For the latter, care must be taken to maintain the application of the VELO_NORMAL condition after a remesh. This condition is applied only to one side of the interface and depends on the ss_to_blks connectivity structure; it may be necessary to force its application, especially after remeshes. To be sure that the proper
  set of conditions is being applied, look at the BC_dup.txt file for nodes along the interface.

--------------
**References**
--------------

GT-003.1: Roll coating templates and tutorial for GOMA and SEAMS, February 29,
2000, P. R. Schunk and Matt Stay

GT-006.3: Slot and Roll coating with remeshing templates and tutorial for GOMA and CUBIT/MAPVAR, August 3, 1999, R. R. Lober and P. R. Schunk

..
	TODO - Line 55 has an image that needs to be replaced with the correct equation.