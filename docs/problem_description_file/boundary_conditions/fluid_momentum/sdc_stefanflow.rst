******************
**SDC_STEFANFLOW**
******************

::

	BC = SDC_STEFANFLOW SS <bc_id> <integer> {char_string}

-----------------------
**Description / Usage**
-----------------------

**(SIC/MOMENTUM)**

This boundary condition represents the specification of the normal component of the
interfacial velocity on one side of the interface. These are *DVI_SIDTIE_VD* boundary
conditions (Moffat, 2001) that have an additional property. The first time encountered
in the formation of the residual, the results of a sub calculation are stored either at the
node structure level or at the surface Gauss point level. The surface reaction and
surface species are specified as part of a surface domain within Chemkin.

The *SURFDOMAINCHEMKIN_STEFAN_FLOW* (shortened to *SDC_STEFANFLOW*
in the *name2* member of the *BC_descriptions* struct in mm_names.h) boundary
condition solves the following equation representing Stefan flow at a boundary.

.. figure:: /figures/113_goma_physics.png
	:align: center
	:width: 90%

where :math:`n_l` is the outward facing normal to the liquid material, :math:`p^l` is the liquid density, :math:`u^l`
is the (mass average) velocity at the current surface quadrature point, and 
:math:`u_s` the
velocity of the mesh (i.e., the interface if the mesh is fixed at the interface). The
summation over *N* species is for the product of molecular weight ( :math:`W_k` ) and the source
term for creation of species k in the liquid ( :math:`S_k^l` ). Note, while it may seem that one side
of the interface is getting special treatment, the combination of this boundary condition
with the KINEMATIC_CHEM boundary condition actually creates a symmetric treatment
of the boundary condition. *SDC_STEFANFLOW* is linked to the SDC_SPECIES_RXN
boundary conditions just as the KINEMATIC_CHEM boundary conditions are by the
expression for the interface reaction. The sum is over all of the interfacial source terms
for species in the phase.

Definitions of the input parameters are as follows:

+------------------+------------------------------------------------------------+
|**SDC_STEFANFLOW**| Name of the boundary condition (<bc_name>).                |
+------------------+------------------------------------------------------------+
|**SS**            | Type of boundary condition (<bc_type>), where **SS**       |
|                  | denotes side set in the EXODUS II database.                |
+------------------+------------------------------------------------------------+
|<bc_id>           | The boundary flag identifier, an integer associated with   |
|                  | <bc_type> that identifies the boundary location (side set  |
|                  | in EXODUS II) in the problem domain.                       |
+------------------+------------------------------------------------------------+
|<integer>         | Element Block ID of the phase on whose side of the         |
|                  | interface this boundary condition will be applied.         |
+------------------+------------------------------------------------------------+
|char_string       | :math:`S_k^l`, string indicating where the surface         |
|                  | source term information for this boundary condition will be|
|                  | obtained. Three options exist:                             |
|                  |                                                            |
|                  |   * **IS_EQUIL_PSEUDORXN**                                 |
|                  |   * **VL_EQUIL_PSEUDORXN**                                 |
|                  |   * **SDC_SURFRXN**                                        |
|                  |                                                            |
|                  | These are boundary conditions that apply to the *Species   |
|                  | Equations*. The last boundary condition is not yet         |
|                  | implemented, so **SDC_SURFRXN** currently does nothing.    |
+------------------+------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample input card:
::

    BC = SDC_STEFANFLOW SS 1   0 VL_EQUIL_PSEUDORXN

The above card will create a strongly integrated boundary condition specifying the
normal component of the velocity on side set 1 on the element block 0 side of the
interface. The source term to be used will be taken from multiple previously specified
*VL_EQUIL_PSEUDORXN* cards.

-------------------------
**Technical Discussion**
-------------------------

* Currently, this card has only been tested on internal interfaces containing
  discontinuous interfaces using the *VL_EQUIL_PSEUDORXN* source term. The
  *SDC_SURFRXN* boundary condition has not been implemented yet.

* The *DVI_SIDTIE_VD* variable is a nomenclature adopted by Moffat (2001) in his
  development of a revised discontinuous variable implementation for *Goma*. It
  pertains to Discontinuous Variable Interfaces (**DVI**) and the strongly integrated
  Dirichlet (**SID**) boundary conditions prescribing the discontinuous value of
  variables on either side of an interface (**TIE** boundary conditions). The user is
  referred to Moffat (2001) for detailed presentation on discontinuous variables.



--------------
**References**
--------------

GTM-015.1: Implementation Plan for Upgrading Boundary Conditions at
Discontinuous-Variable Interfaces, January 8, 2001, H. K. Moffat

.. TODO - Line 26 have photos that needs to be replaced with the real equation.