****************************
**Species Time Integration**
****************************

::

   Species Time Integration= {model_name} <species>

-----------------------
**Description / Usage**
-----------------------

Sharp gradients are often a feature of convective-diffusive computations involving
species. Traditional Galerkin time integration is not optimal under these circumstances.
This optional card is used to change the species time integration scheme to be different
from the global time integration. Each species equation can use a different time
integration. The new time integration schemes are based upon a Taylor-Galerkin
formulation which has better behavior when sharp fronts are present.

Following are the {model_name} options for species time integration, each of which
requires only a species designation to which the model should be applied:

+-----------------------+-------------------------------------------------------------------------------------+
|**STANDARD**           |the input deck formulation, i.e., the global time integration scheme; this is the    |
|                       |default.                                                                             |
|                       |                                                                                     |
|                       | * <species> - the index of the species equation.                                    |
+-----------------------+-------------------------------------------------------------------------------------+
|**TAYLOR_GALERKIN**    |an implicit or semi-implicit Taylor-Galerkin time integration scheme                 |
|                       |                                                                                     |
|                       | * <species> - the index of the species equation.                                    |
+-----------------------+-------------------------------------------------------------------------------------+
|**TAYLOR_GALERKIN_EXP**|An explicit Taylor-Galerkin time integration scheme                                  |
|                       |                                                                                     |
|                       | * <species> - the index of the species equation.                                    |
+-----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

The following sample input card invokes the explicit Taylor-Galerkin time integration of the species equation.

::

   Species Time Integration = TAYLOR_GALERKIN_EXP 0

-------------------------
**Technical Discussion**
-------------------------

The Taylor-Galerkin schemes are designed for advection dominated problems with
sharp fronts where rigorous mass conservation is important.

* **TAYLOR_GALERKIN** uses an implicit or semi-implicit form of the Taylor-
  Galerkin time integrals depending on what is chosen in the input deck.

* **TAYLOR_GALERKIN_EXP** uses an explicit form of the equations and is
  favored for volume-of-fluid simulations where the diffusive character of the
  implicit solver creates mass balance errors. The drawback of explicit time
  integration methods is that the time step used is governed by the Courant limit and
  must be quite small for stability.



--------------
**References**
--------------

No References.