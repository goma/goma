**************************
Pressure Stabilization
**************************

::

	Pressure Stabilization = {yes | no | local | pspp | pspp_e}

-----------------------
Description / Usage
-----------------------

This optional card indicates whether or not pressure stabilization should be used. Valid
options are


yes
    Use the Galerkin Least square pressure stabilization method developed by
    Hughes, et. al. (1986).
local
    Use the Galerkin Least square pressure stabilization method with local
    scaling.
pspp
    Use polynomial stabilized pressure projection stabilization method
    developed by Dohrmann and Bochev (2004). Please see Level Set PSPP
    filtering card if using with the level-set front tracking technique.
pspp_e
    Use polynomial stabilized pressure projection method with upgrade for
    nonuniform/graded meshes (recommended)
no
    Do not use any pressure stabilization.


The amount of pressure stabilization to use is specified with the *Pressure Stabilization
Scaling* card.

The default is **no**, to not use pressure stabilization.

------------
Examples
------------

Following is a sample card:
::

	Pressure Stabilization = yes

-------------------------
Technical Discussion
-------------------------

If input for this card is **yes**, the Hughes, et al. (1986) method adds the residual of the
momentum equation weighted by the gradient of the Galerkin weight function to the
Galerkin continuity equation. The result is that the continuity equation now has a
diagonal term to stabilize it and improve the condition of the matrix, allowing for the
use of iterative solvers. When pressure stabilization is used, equal-order interpolation
can (and should) be used for velocity and pressure, e.g., velocity and pressure both Q2
or both Q1. If input for this card is **no**, then the standard Galerkin finite-element weight
functions are used and velocity and pressure interpolations should be chosen to satisfy
the Babuska-Brezzi condition, e.g., velocity Q2 and pressure Q1 or P1, or velocity Q1
and pressure P0.

An improvement on the Hughes approach was developed by Bochev and Dohrmann
(2004) called the polynomial stabilized pressure projection. In its fundamental form, it
is like PSPG just an additional term on the continuity equation residual that helps
stabilize the pressure, and it is predicated on the fact that the pressure field is governed
by an elliptical equation known as the pressure Poisson equation. Please consult this
paper for details. An additional improvement to that technique was developed
internally to Sandia which better accommodates graded meshes. This technique is
invoked with the pspp_e option, which we recommend.

--------------
References
--------------

Hughes, T. J. R., L. P. Franca and M. Balestra, “A New Finite Element Formulation for
Computational Fluid Dynamics: V. Circumventing the Babuska-Brezzi Condition: A
Stable Petrov-Galerkin Formulation of the Stokes Problem Accommodating Equal-
Order Interpolations,” *Comput. Methods Appl. Mech. Engrg.*, 59 (1986) 85-99.
