Source Terms
################

Source term models cover the internal generation of pressure (in fluids and solids), energy,
species component mass and electrical potential for each of the main differential equations.
Several representations are available for fluids, and the user should be aware of some
consistencies required with density models (see details below). For all of the source models, the
user must heed the following warning:

**Make sure the equation term multipliers for the source terms being used are set to unity**
*(Section 4.12 - Problem Description and Equation specification in Volume 1).*

.. include:: source_terms/navier_stokes_source.rst

.. include:: source_terms/solid_body_source.rst

.. include:: source_terms/mass_source.rst

.. include:: source_terms/heat_source.rst

.. include:: source_terms/species_source.rst

.. include:: source_terms/current_source.rst

.. include:: source_terms/moment_source.rst

.. include:: source_terms/initialize.rst

