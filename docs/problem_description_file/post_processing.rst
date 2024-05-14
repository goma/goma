Post Processing Specifications
##################################

This section lists the post-processing options that are accessible within *Goma*. Each card below triggers calculations of the nodal values of a given function, which are then written to the EXODUS II output file. Normally these values are smoothed before writing them to the output file. For most of these cards a keyword is the only input; if the keyword is **yes**, the post-processing variable is calculated and written to the file; if the keyword is **no**, no output is generated for that variable. All of these cards are optional and can appear in any order.

The sections below list the post-processing options and a brief description of each. *Users are cautioned - for large, time-dependent runs, the output of many post-processing variables may lead to excessively large EXODUS II output files*.

.. toctree::
   :maxdepth: 1

   post_processing/stream_function
   post_processing/streamwise_normal_stress
   post_processing/cross_stream_shear_rate
   post_processing/mean_shear_rate
   post_processing/pressure_contours
   post_processing/fill_contours
   post_processing/concentration_contours
   post_processing/stress_contours
   post_processing/first_invariant_of_strain
   post_processing/second_invariant_of_strain
   post_processing/third_invariant_of_strain
   post_processing/velocity_divergence
   post_processing/particle_velocity_divergence
   post_processing/total_velocity_divergence
   post_processing/electric_field
   post_processing/electric_field_magnitude
   post_processing/enormsq_field
   post_processing/enormsq_field_norm
   post_processing/viscosity
   post_processing/density
   post_processing/lame_mu
   post_processing/lame_lambda
   post_processing/von_mises_strain
   post_processing/von_mises_stress
   post_processing/navier_stokes_residuals
   post_processing/moving_mesh_residuals
   post_processing/mass_diffusion_vectors
   post_processing/diffusive_mass_flux_vectors
   post_processing/mass_fluxlines
   post_processing/energy_conduction_vectors
   post_processing/energy_fluxlines
   post_processing/time_derivatives
   post_processing/mesh_stress_tensor
   post_processing/real_solid_stress_tensor
   post_processing/mesh_strain_tensor
   post_processing/viscoplastic_def_grad_tensor
   post_processing/lagrangian_convection
   post_processing/normal_and_tangent_vectors
   post_processing/error_zz_velocity
   post_processing/error_zz_heat_flux
   post_processing/error_zz_pressure
   post_processing/user_defined_post_processing
   post_processing/porous_saturation
   post_processing/total_density_of_solvents_in_porous_media
   post_processing/density_of_solvents_in_gas_phase_in_porous_media
   post_processing/density_of_liquid_phase_in_porous_media
   post_processing/gas_phase_darcy_velocity_in_porous_media
   post_processing/liquid_phase_darcy_velocity_in_porous_media
   post_processing/capillary_pressure_in_porous_media
   post_processing/grid_peclet_number_in_porous_media
   post_processing/supg_velocity_in_porous_media
   post_processing/vorticity_vector
   post_processing/map_conf_stress
   post_processing/viscous_stress
   post_processing/fluid_stress
