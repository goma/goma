********************
**SHELL_LUBP_SOLID**
********************

::

	BC = SH_LUBP_SOLID SS <bc_id> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/R_MESH1/R_MESH2/RMESH3)**

This vector boundary condition card balances the stress in an abutting continuum elastic solid with the lubrication forces (pressure and shear) in a surface shell. The boundary condition is applied to a sideset. Please see notes below on the sideset features which must be specified.

============= ===========================================================
SH_LUBP_SOLID Name of boundary condition.
SS            Type of boundary condition (<bc_type>), where
              SS denotes sideset in the EXODUS II database.
<bc_id>       The boundary flag identifier, an integer associated with
              <bc_type> that identifies the boundary location (nodeset
              in EXODUS II) in the problem domain.
<float1>      Scaling factor. Normally set this to 1.0, unless a 
              stressbalance scale is required due to 
              nondimensionalization.
============= ===========================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = SH_LUBP_SOLID SS   100 1.0

This boundary condition is applied at sideset 100.

-------------------------
**Technical Discussion**
-------------------------

* The mathematical form of the boundary condition is

.. math::

   \boldsymbol{n}\cdot\boldsymbol{\sigma}=\boldsymbol{k}P_{lub}+\boldsymbol{i}\left(\frac{h}{12}\frac{\partial P}{\partial x}+\frac{u\mu}{h}\right)+\boldsymbol{j}\left(\frac{h}{12}\frac{\partial P}{\partial x}+\frac{v\mu}{h}\right)

* This condition is similar to FLUID_SOLID and SOLID_FLUID boundary
  conditions for the case of fluid-structure interaction between two continuum
  regions, one fluid and one solid.

* Note that the sideset as generated in CUBIT or related software is actually 
  attached to the continuum domain and not the shell face, as those faces (top and bottom of sheet and not the edges) are not true finite element sides. Most mesh generators will not allow sidesets to be include shell element faces. GOMA figures out the right thing to do.
