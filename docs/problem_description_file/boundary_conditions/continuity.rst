~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category 7: Continuity Equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The continuity equation rarely requires a boundary condition as it represents an overall mass
balance constraint on the velocity field for the fluid, viz. normally it is used to enforce
incompressibility. Boundary conditions for pressure are most often put on the fluid-momentum
equations as a part of the stress condition at an inflow or outflow plane (see for example boundary
condition cards *FLOW_PRESSURE, FLOW_HYDROSTATIC*, etc. ). On occasion, however, we
can use a pressure condition as a pressure datum, as the Dirichlet pressure condition below
allows, though the user must keep in mind that it is a condition on continuity and not momentum.
When using pressure stabilization, viz. PSPG techniques, then also there is an occasional need for
a boundary condition on this equation.

.. include:: /problem_description_file/boundary_conditions/continuity/p.rst

.. include:: /problem_description_file/boundary_conditions/continuity/pspg.rst

.. include:: /problem_description_file/boundary_conditions/continuity/pressure_datum.rst

