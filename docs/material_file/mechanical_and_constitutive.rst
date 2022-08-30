Mechanical Properties and Constitutive Equations
####################################################

This section of the material property input specifies the type of model, for both solids and fluids,
that relates stress and strain (or strain-rate) as well as the various parameters for these models.
Models for solids are relatively simple compared to solid mechanics codes but cover the primary
needs in fluid-solid problems. The models for fluids are quite extensive, covering Newtonian,
generalized-Newtonian, rate-dependent models, thermally-dependent models, curing and particleladen
models and combinations of these. These properties are used in the solid and fluid
momentum conservation equations.

.. toctree::
   :maxdepth: 1

   mechanical_and_constitutive/solid_constitutive_equation
   mechanical_and_constitutive/plasticity_equation
   mechanical_and_constitutive/convective_lagrangian_velocity
   mechanical_and_constitutive/lame_mu
   mechanical_and_constitutive/lame_lambda
   mechanical_and_constitutive/stress_free_solvent_vol_frac
   mechanical_and_constitutive/solid_thermal_expansion
   mechanical_and_constitutive/solid_reference_temperature
   mechanical_and_constitutive/plastic_viscosity
   mechanical_and_constitutive/evp_yield_stress
   mechanical_and_constitutive/pseudo_solid_constitutive_equation
   mechanical_and_constitutive/pseudo_solid_lame_mu
   mechanical_and_constitutive/pseudo_solid_lame_lambda
   mechanical_and_constitutive/liquid_constitutive_equation
   mechanical_and_constitutive/viscosity
   mechanical_and_constitutive/low_rate_viscosity
   mechanical_and_constitutive/power_law_exponent
   mechanical_and_constitutive/high_rate_viscosity
   mechanical_and_constitutive/time_constant
   mechanical_and_constitutive/aexp
   mechanical_and_constitutive/thermal_exponent
   mechanical_and_constitutive/thermal_wlf_constant2
   mechanical_and_constitutive/yield_stress
   mechanical_and_constitutive/yield_exponent
   mechanical_and_constitutive/suspension_maximum_packing
   mechanical_and_constitutive/suspension_species_number
   mechanical_and_constitutive/cure_gel_point
   mechanical_and_constitutive/cure_a_exponent
   mechanical_and_constitutive/cure_b_exponent
   mechanical_and_constitutive/cure_species_number
   mechanical_and_constitutive/unreacted_gel_temperature
   mechanical_and_constitutive/polymer_constitutive_equation
   mechanical_and_constitutive/ptt_form
   mechanical_and_constitutive/polymer_stress_formulation
   mechanical_and_constitutive/polymer_weight_function
   mechanical_and_constitutive/polymer_shift_function
   mechanical_and_constitutive/polymer_weighting
   mechanical_and_constitutive/polymer_shock_capturing
   mechanical_and_constitutive/discontinuous_jacobian_formulation
   mechanical_and_constitutive/adaptive_viscosity_scaling
   mechanical_and_constitutive/polymer_viscosity
   mechanical_and_constitutive/polymer_time_constant
   mechanical_and_constitutive/polymer_yield_stress
   mechanical_and_constitutive/mobility_parameter
   mechanical_and_constitutive/ptt_xi_parameter
   mechanical_and_constitutive/ptt_epsilon_parameter
   mechanical_and_constitutive/surface_tension
   mechanical_and_constitutive/second_level_set_conductivity
   mechanical_and_constitutive/second_level_set_density
   mechanical_and_constitutive/second_level_set_heat_capacity
   mechanical_and_constitutive/second_level_set_momentum_source
   mechanical_and_constitutive/second_level_set_viscosity
   mechanical_and_constitutive/shell_bending_stiffness
