****************
LS_CAP_HYSING
****************

::

	BC = LS_CAP_HYSING LS <double>

-----------------------
Description / Usage
-----------------------

**(EMB/VECTOR MOMENTUM)**

Embedded boundary condition for applying a surface tension source term
with level set.  Has the added advantage of a time scaled diffusion
term to aid with reducing spurious currents.The surface tension value
used in this boundary condition is obtained from the Surface Tension
material parameter defined in the mat file.

Not compatible with subelement integration.

A description of the input parameters follows:

+--------------------+---------------------------------------------+
|LS_CAP_HYSING       |Boundary condition name                      |
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

   BC = LS_CAP_HYSING LS 1.0

-------------------------
Technical Discussion
-------------------------

Adds the following embedded boundary condition to the momentum equations:

.. math::

   - \int \nabla \mathbf{v} : (\mathbf{I} - \mathbf{\hat{n}}^n \mathbf{\hat{n}}^n) d\mathbf{x} 
   - \beta \Delta t^{n+1} \int \nabla_s \mathbf{v} : (\sigma \delta_{\alpha}(\phi) \nabla_s \mathbf{u}^{n+1}) d\mathbf{x}

Where, :math:`\beta` is a scaling term, :math:`n` represents a timestep, and :math:`\mathbf{\hat{n}}` represents the normal
found from the level set function, and :math:`\mathbf{u}` is the velocity.

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

Hysing, S.R., Turek, S., Kuzmin, D., Parolini, N., Burman, E.,
Ganesan, S. and Tobiska, L., 2009. Quantitative benchmark computations
of two‐dimensional bubble dynamics. International Journal for
Numerical Methods in Fluids, 60(11), pp.1259-1288.
