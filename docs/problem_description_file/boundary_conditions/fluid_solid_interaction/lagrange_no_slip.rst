********************
**LAGRANGE_NO_SLIP**
********************

::

	BC = LAGRANGE_NO_SLIP SS <bc_id> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(CONTACT_SURF/R_LAGR_MULT1)**

This boundary condition is used to apply a kinematic Lagrange multiplier constraint to
a solid/fluid boundary while using Goma’s overset grid capability. The condition is
used when the complete fluid-structure interaction problem is being solved, viz.
stresses between fluid and solid are both accommodated as is the dynamics of the
structure and fluid. In contrast, *Goma* allows for a structure to be moved through a
fluid under prescribed kinematics, and in that case a different Lagrange multiplier
constraint is advocated (see LS_NO_SLIP, for example). Two integer inputs together
with a sideset ID integer are required for this boundary condition:

===================== ==============================================================
**LAGRANGE_NO_SLIP**  Name of the boundary condition.
**SS**                Type of boundary condition (<bc_type>), where **NS** denotes
                      node set in the EXODUS II database.
<bc_id>               The boundary flag identifier, an integer associated with
                      <bc_type> that identifies the boundary location (node set in
                      EXODUS II) in the problem domain.
<integer1>            Element block ID of solid phase from the EXODUS II
                      database.
<integer2>            Element block ID of liquid phase from the EXODUS II
                      database.
===================== ==============================================================

------------
**Examples**
------------

Following is a sample card:
::

   BC = LAGRANGE_NO_SLIP SS 2 1 2

In this case the kinematic condition (viz. a velocity match of fluid and solid at the
interface) is applied to the interface imprinted by sideset 2. That side set is imprinted
on the background fluid mesh. The solid material block ID is 1 in this case and the
background fluid material ID is 2. 

-------------------------
**Technical Discussion**
-------------------------

In this work, the governing equations consist of a fluid momentum balance:

.. figure:: /figures/188_goma_physics.png
	:align: center
	:width: 90%

a mass balance:

.. figure:: /figures/189_goma_physics.png
	:align: center
	:width: 90%

and a solid momentum balance:

.. figure:: /figures/190_goma_physics.png
	:align: center
	:width: 90%

The kinematic constraint at the fluid-solid interface is:

.. figure:: /figures/191_goma_physics.png
	:align: center
	:width: 90%

and the level set function is evaluated at each fluid mesh node by:

.. figure:: /figures/192_goma_physics.png
	:align: center
	:width: 90%

The first four equations are written in a Galerkin/Finite form, with :math:`\phi_i` representing the
weighting functions at node *i*. The first three equations are enforced at all nodes *i* that
contain the appropriate degrees of freedom (viz. solid or fluid dofs). The fourth
equation applies at the solid-liquid interface. :math:`\rho_f` and :math:`\rho_s` are the fluid and solid material
densities, respectively, :math:`\underline{v}` is the fluid velocity, :math:`\underline{F}` represents any body forces such as
gravity, :math:`\tau` is the fluid stress tensor, :math:`\underline{\gamma}` is the Lagrange multiplier vector unknown, :math:`\underline{\chi}` is the
solid displacement vector unknown, :math:`\sigma` is the solid stress tensor, *f* is the level set
unknown, :math:`\Theta` is a step function which is -1 for points within the region occupied by the
solid and +1 outside this region, :math:`\underline{x_i}` and :math:`\underline{x_s}` are the position vectors of a fluid node and of
the closest point to it on the solid boundary, respectively, *V* is the fluid volume domain,
*S* is the solid volume domain, and :math:`\Gamma` is the solid boundary (interface) surface domain.

Noteworthy is that this boundary condition applies the second-to-last “kinematic”
constraint.



--------------
**References**
--------------

GT-026.3: Goma’s Overset Mesh Capability. Randy Schunk.