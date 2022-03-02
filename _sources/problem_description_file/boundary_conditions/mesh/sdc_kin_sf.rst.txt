**************
**SDC_KIN_SF**
**************

::

	BC = SDC_KIN_SF SS <bc_id> <integer> {char_string}

-----------------------
**Description / Usage**
-----------------------

**(SIC/ROTATED MESH)**

This boundary condition represents the specification of the normal component of the
mesh velocity. This is a DVI_MULTI_PHASE_SINGLE boundary condition that has an
additional property. The first time encountered in the formation of the residual, the
results of a subcalculation are stored either at the node structure level or at the surface
gauss point level. The surface reaction and surface species are specified as part of a
surface domain within Chemkin.

The SURFDOMAINCHEMKIN_KIN_STEFAN_FLOW boundary condition (shortened to
SDC_KIN_SF in the *name2* member of the *BC_descriptions* struct in mm_names.h)
solves the following equation representing Stefan flow at a boundary.

.. figure:: /figures/070_goma_physics.png
	:align: center
	:width: 90%

where :math:`n_1` is the outward facing normal to the liquid material, :math:`p^1` is the liquid density, :math:`u^1`
is the (mass average) velocity at the current surface quadrature point, and :math:`u_s` the
velocity of the mesh (i.e., the interface if the mesh is fixed at the interface). The
summation over *N* species is for the product of molecular weight ( :math:`W_k` ) and the source
term for creation of species k in the liquid ( :math:`S^1_k` ). SDC_KIN_SF is linked to the
SDC_SPECIES_RXN boundary conditions just as the KINEMATIC_CHEM boundary
conditions are by the expression for the interface reaction. The sum is over all of the
interfacial source terms for species in the phase.

Definitions of the input parameters are as follows:

+---------------+----------------------------------------------------------------+
|**SDC_KIN_SF** | Name of the boundary condition (<bc_name>).                    |
+---------------+----------------------------------------------------------------+ 
|**SS**         | Type of boundary condition (<bc_type>), where **SS** denotes   |
|               | side set in the EXODUS II database.                            |
+---------------+----------------------------------------------------------------+
|<bc_id>        | The boundary flag identifier, an integer associated with       |
|               | <bc_type> that identifies the boundary location (side set in   |
|               | EXODUS II) in the problem domain.                              |
+---------------+----------------------------------------------------------------+
|<integer>      | lement Block ID of the phase on whose side of the              |
|               | interface this boundary condition will be applied.             |
+---------------+----------------------------------------------------------------+
|char_string    | :math:`S^1_k` string indicating where the surface source term  |
|               | information for this boundary condition will be                |
|               | obtained. Three options exist:                                 |
|               |                                                                |
|               |   * **IS_EQUIL_PSEUDORXN**                                     |
|               |   * **VL_EQUIL_PSEUDORXN**                                     |
|               |   * **SDC_SURFRXN**                                            |
|               |                                                                |
|               | These are boundary conditions that apply to the *Species       |
|               | Equations*. The last boundary condition is not yet             |
|               | implemented, so **SDC_SURFRXN** currently does nothing.        |
+---------------+----------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample card:
::

     BC = SDC_KIN_SF SS 1   0 VL_EQUIL_PSEUDORXN

The above card will create a strongly integrated boundary condition specifying the
normal component of the velocity on side set 1 on the element block 0 side of the
interface. The source term to be used in the above equation will be taken from multiple
previously specified multiple VL_EQUIL_PSEUDORXN cards.

-------------------------
**Technical Discussion**
-------------------------

* This boundary condition is exactly the same as SDC_STEFANFLOW, except for the
  fact that it is applied on the normal component of the mesh velocity instead of the
  normal component of the mass averaged velocity. It is similar to a single phase
  boundary condition, because all of its input comes from one side of the interface.
  Thus, it can equally be applied to external surfaces as well as internal ones with
  some development work.

* Currently, it has only been tested out on internal boundaries using the
  IS_EQUIL_PSEUDORXN source term.

* The DVI_MULTI_PHASE_SINGLE variable is a nomenclature adopted by Moffat
  (2001) in his development of a revised discontinuous variable implementation for
  *Goma*. It pertains to Discontinuous Variable Interfaces (**DVI**) and boundary
  conditions that involve the addition of a surface integral to each side of an internal
  boundary for a variable that is continuous across the interface. The user is referred
  to Moffat (2001) for detailed presentation on discontinuous variables.



--------------
**References**
--------------

GTM-015.1: Implementation Plan for Upgrading Boundary Conditions at
Discontinuous-Variable Interfaces, January 8, 2001, H. K. Moffat

_____________________________________________________________________________

.. 
	TODO - The image in line 26 needs to be replaced with the correct equation.
