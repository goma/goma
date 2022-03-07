~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category 12: Fill Equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The so-called *Fill* equation is used by the volume-of-fluid and level-set Eulerian interface
tracking in *Goma*. Basically it is a statement of Lagrangian invariance and is hence a hyperbolic
statement of the so-called kinematic equation. Given a velocity field, this equation advances the
fill function as a set of material points; hence material surfaces remain ostensibly intact. The
boundary conditions in this section are used to specify the level-of-fill at a boundary at which a
fluid of a specific phase is flowing into the problem domain.

.. include:: /problem_description_file/boundary_conditions/fill/f.rst

.. include:: /problem_description_file/boundary_conditions/fill/fill_inlet.rst
