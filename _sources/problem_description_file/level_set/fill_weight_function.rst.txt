************************
Fill Weight Function
************************

::

	Fill Weight Function = {Galerkin | Taylor-Galerkin | SUPG}

-----------------------
Description / Usage
-----------------------

Sets the weight function used for the *FILL* equation for either the VOF or Level Set
methods. The options for this card are as follows:

Galerkin
    Name of the weight function formulation. This option requests a standard
    Galerkin finite element weighted residual treatment. A floating point
    parameter is not used for this option.

Taylor-Galerkin
    Name of the weight function formulation

SUPG
    Name of the weight function formulation. This option requests a Streamwise
    Upwinding Petrov Galerkin formulation. No floating point parameter is
    required.

The default value for the *Fill Weight Function* is Taylor-Galerkin.

------------
Examples
------------

This is a sample card:
::

	Fill Weight Function = Galerkin

-------------------------
Technical Discussion
-------------------------

This card selects the integration/weight function used in solving for the VOF color
function or the level set distance function (i.e., the *FILL* unknown). The user should
refer to the tutorial on Level Set Computations for a detailed description of level set
interface tracking. (See References.)

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer

A. N. Brooks and T. J. R. Hughes, “Streamline Upwind/Petrov-Galerkin Formulations
for Convection Dominated Flows with Particular Emphasis on the Incompressible
Navier-Stokes Equations,” Comp. Math. In Appl. Mechanics and Eng., 32, 199 - 259
(1992).

A. J. A. Unger, P. A. Forsyth and E. A. Sudicky, “Variable spatial and temporal
weighting schemes for use in multi-phase compositional problems,” Advances in Water
Resources, 19, 1 - 27 (1996).

R. Helmig and R. Huber, “Comparison of Galerkin-type discretization techniques for
two-phase flow in heterogeneous porous media,” Advances in Water Resources, 21,
697-711 (1998).

E. Gundersen and H. P. Langtangen, “Finite Element Methods for Two-Phase Flow in
Heterogeneous Porous Media,” in Numerical Methods and Software Tools in Industrial
Mathematics, Morten Daehlen, Aslak Tveito, Eds., Birkhauser, Boston, 1997.

S. F. Bradford and N. D. Katopodes, “The anti-dissipative, non-monotone behavior of
Petrov-Galerkin Upwinding,” International J. for Numerical Methods in Fluids, v. 33,
583-608 (2000).
