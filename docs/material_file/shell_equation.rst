Shell Equation Properties and Models
########################################

In this section we list all “material-region” specific models and properties associated with
GOMA’s extensive shell equation capability. Currently we have specialized shell equations for
Reynolds lubrication flow (lubp), open Reynolds film flow (shell_film_H), energy (shell_energy,
convection and diffusion, coupled with lubrication), thin porous media (closed cell and open cell),
melting and phase change and more. While many of these cards are actual material properties,
most are geometry and kinematic related. The most appropriate place for these cards are region/
material files because they are actually boundary conditions and related parameters which arise
from the reduction of order (integration through the thin film). For more information, please see
the shell-equation tutorial (GT-036).

.. include:: shell_equation/upper_height_function_constants.rst

.. include:: shell_equation/lower_height_function_constants.rst

.. include:: shell_equation/upper_velocity_function_constants.rst

.. include:: shell_equation/lower_velocity_function_constants.rst

.. include:: shell_equation/upper_contact_angle.rst

.. include:: shell_equation/lower_contact_angle.rst

.. include:: shell_equation/lubrication_fluid_source.rst

.. include:: shell_equation/lubrication_momentum_source.rst

.. include:: shell_equation/turbulent_lubrication_model.rst

.. include:: shell_equation/shell_energy_source_QCONV.rst

.. include:: shell_equation/shell_energy_source_sliding_contact.rst

.. include:: shell_equation/shell_energy_source_viscous_dissipation.rst

.. include:: shell_equation/shell_energy_source_external.rst

.. include:: shell_equation/FSI_deformation_model.rst

.. include:: shell_equation/film_evaporation_model.rst

.. include:: shell_equation/disjoining_pressure_model.rst

.. include:: shell_equation/diffusion_coefficient_model.rst

.. include:: shell_equation/porous_shell_radius.rst

.. include:: shell_equation/porous_shell_height.rst

.. include:: shell_equation/porous_shell_closed_porosity.rst

.. include:: shell_equation/porous_shell_closed_gas_pressure.rst

.. include:: shell_equation/porous_shell_atmospheric_pressure.rst

.. include:: shell_equation/porous_shell_reference_pressure.rst

.. include:: shell_equation/porous_shell_cross_permeability.rst

.. include:: shell_equation/porous_shell_gas_diffusivity.rst

.. include:: shell_equation/porous_shell_gas_temperature_constant.rst

.. include:: shell_equation/porous_shell_henrys_law_constant.rst
