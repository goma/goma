******************
LS_CAP_DENNER_DIFF
******************

::

	BC = LS_CAP_DENNER_DIFF LS <double>

-----------------------
Description / Usage
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition is expected to be paired with another boundary
condition that represents the surface tension term.

e.g.

LS_CAPILLARY
LS_CAP_DIV_N
etc..

And should not be used with LS_CAP_HYSING

Embedded boundary condition for applying a diffusion to the level set
capillary boundary condition.  A time scaled diffusion term to aid
with reducing spurious currents.

Not compatible with subelement integration.

A description of the input parameters follows:

+--------------------+---------------------------------------------+
|LS_CAP_DENNER_DIFF  |Boundary condition name                      |
|                    |                                             |
+--------------------+---------------------------------------------+
|LS                  |Indicator that this is a level set boundary  |
|                    |condition                                    |
|                    |                                             |
+--------------------+---------------------------------------------+
|<double>            |Scaling of the diffusive term                |
|                    |                                             |
+--------------------+---------------------------------------------+

------------
Examples
------------

An example:
::

   BC = LS_CAP_DENNER_DIFF LS 1.0

-------------------------
Technical Discussion
-------------------------

Adds the following embedded boundary condition to the momentum
equations:

.. math::
   - \beta \Delta t^{n+1} \int \nabla_s \mathbf{v} : (\sigma \delta_{\alpha}(\phi) \nabla_s \mathbf{u}^{n+1}) d\mathbf{x}

And      

.. math::
   \nabla_s f = \nabla f - \mathbf{\hat{n}}^n (\mathbf{\hat{n}}^n \cdot \nabla f)

Where, :math:`\beta` is a scaling term, :math:`n` represents a
timestep, and :math:`\mathbf{\hat{n}}` represents the normal found
from the level set function (or from the normal equations when
enabled), and :math:`\mathbf{u}` is the velocity.

----------
Theory
----------

No Theory.

--------
FAQs
--------

No FAQs.

--------------
References
--------------
Denner, F., Evrard, F., Serfaty, R. and van Wachem,
B.G., 2017. Artificial viscosity model to mitigate numerical artefacts
at fluid interfaces with surface tension. Computers & Fluids, 143,
pp.59-72.
