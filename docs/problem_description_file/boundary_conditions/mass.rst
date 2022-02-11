~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category 6: Mass Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The collection of boundary conditions in this category are applied to the mass balance equations,
specifically the species component balance equations. Most boundary conditions are weakly
integrated conditions defining fluxes at internal or external surfaces, although strongly integrated
and Dirichlet conditions are also available to control known values of dependent variables or
integrated quantities. Boundary conditions are available for chemical species as well as charged
species, suspensions and liquid metals. An important capability in *Goma* is represented by the
discontinuous variable boundary conditions, for which users are referred to Schunk and Rao
(1994) and Moffat (2001). Care must be taken if the species concentration is high enough to be
outside of the *dilute species* assumption, in which case transport of species through boundaries
will affect the volume of the bounding fluids. In these cases, users are referred to the
*VNORM_LEAK* condition for the fluid momentum equations and to *KIN_LEAK* for the solid
momentum (mesh) equations. And finally, users are cautioned about different bases for
concentration (volume, mass, molar) and several discussions on or references to units.

.. include:: /problem_description_file/boundary_conditions/mass/y.rst

.. include:: /problem_description_file/boundary_conditions/mass/yuser.rst

.. include:: /problem_description_file/boundary_conditions/mass/y_discontinuous.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_const.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_equil.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_sulfidation.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_sus.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_bv.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_hor.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_orr.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_user.rst

.. include:: /problem_description_file/boundary_conditions/mass/yflux_alloy.rst

.. include:: /problem_description_file/boundary_conditions/mass/ytotalflux_const.rst

.. include:: /problem_description_file/boundary_conditions/mass/vl_equil.rst

.. include:: /problem_description_file/boundary_conditions/mass/vl_poly.rst

.. include:: /problem_description_file/boundary_conditions/mass/vl_equil_pseudorxn.rst

.. include:: /problem_description_file/boundary_conditions/mass/is_equil_pseudorxn.rst

.. include:: /problem_description_file/boundary_conditions/mass/surface_charge.rst
