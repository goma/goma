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

.. toctree::
   :maxdepth: 1

   shell_equation/upper_height_function_constants
   shell_equation/lower_height_function_constants
   shell_equation/upper_velocity_function_constants
   shell_equation/lower_velocity_function_constants
   shell_equation/upper_contact_angle
   shell_equation/lower_contact_angle
   shell_equation/lubrication_fluid_source
   shell_equation/lubrication_momentum_source
   shell_equation/turbulent_lubrication_model
   shell_equation/shell_energy_source_QCONV
   shell_equation/shell_energy_source_sliding_contact
   shell_equation/shell_energy_source_viscous_dissipation
   shell_equation/shell_energy_source_external
   shell_equation/FSI_deformation_model
   shell_equation/film_evaporation_model
   shell_equation/disjoining_pressure_model
   shell_equation/diffusion_coefficient_model
   shell_equation/porous_shell_radius
   shell_equation/porous_shell_height
   shell_equation/porous_shell_closed_porosity
   shell_equation/porous_shell_closed_gas_pressure
   shell_equation/porous_shell_atmospheric_pressure
   shell_equation/porous_shell_reference_pressure
   shell_equation/porous_shell_cross_permeability
   shell_equation/porous_shell_gas_diffusivity
   shell_equation/porous_shell_gas_temperature_constant
   shell_equation/porous_shell_henrys_law_constant
