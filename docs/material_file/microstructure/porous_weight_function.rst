**************************
**Porous Weight Function**
**************************

::

   Porous Weight Function = {GALERKIN | SUPG} <float>

-----------------------
**Description / Usage**
-----------------------

This required card is used to specify the weight function form on the capacitance term
of the Darcy flow equations for partially saturated flow (viz. for Media Type
specifications of **POROUS_PART_SAT** and **POROUS_UNSAT**, and
**POROUS_TWO_PHASE**.) The standard approach is to use a Galerkin formulation,
but often times the SUPG option allows for a more stable time integration algorithm
using the classic Streamwise Upwinding Petrov Galerkin weight function (see
references below). The model options for this card are as follows:

+-------------------+-------------------------------------------------------------------------------------+
|**GALERKIN**       |Name of the weight function formulation. This option requests a standard Galerkin    |
|                   |finite element weighted residual treatment. A parameter is required, viz. <float>,   |
|                   |but it is not used by *Goma*; it should be set to zero.                              |
|                   |                                                                                     |
|                   | * <float> - 0.0                                                                     |
+-------------------+-------------------------------------------------------------------------------------+
|**SUGP**           |Name of the weight function formulation. This option requests a streamwise upwinding |
|                   |Petrov-Galerkin formulation. A floating point parameter is required as a             |
|                   |SUPG weighting parameter and it should be set between 0.0 (for no upwinding) and 1.0 |
|                   |(for full upwinding).                                                                |
|                   |                                                                                     |
|                   | * <float> - a SUPG weighting parameter                                              |
+-------------------+-------------------------------------------------------------------------------------+

The default model if this card is missing is **GALERKIN**.

------------
**Examples**
------------

An example card

::

   Porous Weight Function = SUPG 1.0

-------------------------
**Technical Discussion**
-------------------------

As mentioned above, this card is used to invoke a streamwise upwinding scheme for
purposes of stabilizing the solution around steep saturation fronts. Galerkin finite
element treatment is often an extremely inaccurate discretization for propagating a
discontinuity, such as is the case around these fronts, and often has to be supplemented
with streamwise diffusion and/or mass lumping so that the saturation variable remains
monotonic and well behaved, viz. to keep it from going below zero. Another expedient
to aid in keeping the front smooth and monotonic is to use mass lumping (cf. *Mass
Lumping* card).



--------------
**References**
--------------

GTM-029.0: SUPG Formulation for the Porous Flow Equations in Goma, H. K.
Moffat, August 2001 (DRAFT).

Bradford, S. F. and N. D. Katopodes, “The anti-dissipative, non-monotone behavior of
Petrov-Galerkin Upwinding,” International J. for Numerical Methods in Fluids, v. 33,
583-608 (2000).

Brooks, A. N. and T. J. R. Hughes, “Streamline Upwind/Petrov-Galerkin Formulations
for Convection Dominated Flows with Particular Emphasis on the Incompressible
Navier-Stokes Equations,” Comp. Math. In Appl. Mechanics and Eng., 32, 199 - 259
(1992).

Gundersen, E. and H. P. Langtangen, “Finite Element Methods for Two-Phase Flow in
Heterogeneous Porous Media,” in Numerical Methods and Software Tools in Industrial
Mathematics, Morten Daehlen, Aslak Tveito, Eds., Birkhauser, Boston, 1997.

Helmig, R. and R. Huber, “Comparison of Galerkin-type discretization techniques for
two-phase flow in heterogeneous porous media,” Advances in Water Resources, 21,
697-711 (1998).

Unger, A. J. A., P. A. Forsyth and E. A. Sudicky, “Variable spatial and temporal
weighting schemes for use in multi-phase compositional problems,” Advances in Water
Resources, 19, 1 - 27 (1996).