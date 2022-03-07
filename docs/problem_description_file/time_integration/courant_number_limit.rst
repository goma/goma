************************
Courant Number Limit
************************

::

	Courant Number Limit = <float>

-----------------------
Description / Usage
-----------------------

This parameter’s roll is to control time step growth based on the well-known Courant
number criterion. This card applies only to level-set problems. This card imposes an
upper limit on the time step size, irrespective of the variable time integrator already in
place.

<float>
    Any floating point number to indicate the Courant number limit.

------------
Examples
------------

A sample card that sets the Courant number to 0.2 is:
::

	Courant Number Limit = 0.2

-------------------------
Technical Discussion
-------------------------

See GT-034 for a thorough discussion.

----------
Theory
----------

The time step limit imposed by this limit is computed as

.. math:: 

   \mathrm{d}t_{\mathrm{limit}} = C \min_{e} \left| \frac{h_e}{ \lVert \hat{U} \rVert_e} \right|

Here :math:`e` is the element, :math:`h_e` is the average size of the element, :math:`C` is the specified Courant
number, and

.. math:: 

   \lVert \hat{U} \rVert_e = \frac{\int_{e}^{} \delta_{\alpha} \left( \phi \right) \underline{v} \cdot \underline{n} \, d \Omega}{\int_{e}^{} \delta_{\alpha} \left( \phi \right) \, d \Omega}


--------------
References
--------------

GT-034: Tutorial on time step parameter selection for level-set problems in GOMA.
April 1, 2006. D. R. Noble
