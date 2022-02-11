**************************************
**Discontinuous Jacobian Formulation**
**************************************

::

   Discontinuous Jacobian Formulation = {model_name} <float>

-----------------------
**Description / Usage**
-----------------------

This optional card is used to specify the off element Jacobian contributions for the
discontinuous Galerkin (DG) discretization of the polymer stress equation. These terms
are important because the DG method uses stress information from upstream elements
to determine the flux in the current element. If the off element Jacobians are not
included, convergence is poor, but including these terms greatly increases the
complexity of the code, the matrix bandwidth and the matrix solution time.

The default sets this option to false, implying that no off element Jacobians are
included. Valid options for {model_name} are:

+-----------------+------------------------------------------------------------------------------------------------------------+
|**FULL**         |adds in the full complement of off-element Jacobians; no floating point data required. This option does not |
|                 |always work in parallel computations.                                                                       |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**EXPLICIT**     |approximates the off-element Jacobians by adding terms to the residual equation based on the previous       |
|                 |iteration.                                                                                                  |
|                 |                                                                                                            |
|                 | * <float> - scales the lagged term.                                                                        |
+-----------------+------------------------------------------------------------------------------------------------------------+
|**SEGREGATED**   |approximates the off-element Jacobians by adding terms to the residual equation based on a mass lumping at  |
|                 |the current iteration.                                                                                      |
|                 |                                                                                                            |
|                 | * <float> - scales the lumped term.                                                                        |
+-----------------+------------------------------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that set the discontinuous Jacobian formulation to full.

::

   Discontinuous Jacobian Formulation = FULL

The following is a sample card that set the discontinuous Jacobian formulation to
explicit. Note this is more of a research option than a production one and the choice of
scaling requires tuning for each problem.

::

   Discontinuous Jacobian Formulation = EXPLICIT 0.1

The following is a sample card that set the discontinuous Jacobian formulation to
segregated. Note this is more of a research option than a production one and the choice
of scaling requires tuning for each problem.

::

   Discontinuous Jacobian Formulation = SEGREGATED 0.2

-------------------------
**Technical Discussion**
-------------------------

For a discussion of the discontinuous Galerkin method see Fortin and Fortin (1989),
Baaijens (1994) or Baaijens (1998). Internal (Sandia) users may find T. A. Baer’s
Gordon Conference presentation (1997) helpful.



--------------
**References**
--------------

Baaijens, F. P. T. , “Application of Low-Order Discontinuous Galerkin Method to the
Analysis of Viscoelastic Flows,” J. Non-Newtonian Fluid Mech., 52, 37-57 (1994).

Baaijens, F. P. T., “An Iterative Solver for the DEVSS/DG Method with Application to
Smooth and Non-smooth Flows of the Upper Convected Maxwell Fluid,” J. Non-
Newtonian Fluid Mech., 75, 119-138 (1998).

Fortin, M. and A. Fortin, “A New Approach for the FEM Simulations of Viscoelastic
Flow, J. Non-Newtonian Fluid Mech., 32, 295-310 (1989).


