~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category 8: Porous Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following conditions are applied as boundary conditions to the porous-flow equations. These
conditions include strong Dirichlet conditions, such as hard sets on porous phase pressure on a
boundary as a constant or function of position, weak-form conditions, such as a specified phase
flux from a convective mass transfer model or a constant flux, and a host of interfacial conditions
for impregnation, etc. The porous flow equations are actually scalar equations that represent
component mass balances. Specifically, there is one component mass balance for the liquid phase,
one for the gas phase, and one for the solid phase. The corresponding three dependent variables in
these balances are the liquid phase pressure, the gas phase pressure, and the porosity, respectively.
These variables are related to the flow through a boundary by their normal gradients (Darcy’s law
formulation) and to the local inventory of liquid and gas through the saturation function. These
implicit terms can often lead to some confusion in setting the boundary conditions so it is
recommended that the user consult the supplementary documentation referenced in the following
porous boundary condition cards.

.. include:: /problem_description_file/boundary_conditions/porous/porous_liq_pressure.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_liq_flux_const.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_gas_pressure.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_gas.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_gas_flux_const.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_conv.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_flux.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_pressure.rst

.. include:: /problem_description_file/boundary_conditions/porous/p_liq_user.rst

.. include:: /problem_description_file/boundary_conditions/porous/porous_temperature.rst
