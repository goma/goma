Post Processing Specifications
##################################

This section lists the post-processing options that are accessible within *Goma*. Each card below triggers calculations of the nodal values of a given function, which are then written to the EXODUS II output file. Normally these values are smoothed before writing them to the output file. For most of these cards a keyword is the only input; if the keyword is **yes**, the post-processing variable is calculated and written to the file; if the keyword is **no**, no output is generated for that variable. All of these cards are optional and can appear in any order.

The sections below list the post-processing options and a brief description of each. *Users are cautioned - for large, time-dependent runs, the output of many post-processing variables may lead to excessively large EXODUS II output files*.

.. include:: post_processing/stream_function.rst

.. include:: post_processing/streamwise_normal_stress.rst

.. include:: post_processing/cross_stream_shear_rate.rst

.. include:: post_processing/mean_shear_rate.rst

.. include:: post_processing/pressure_contours.rst

.. include:: post_processing/fill_contours.rst

.. include:: post_processing/concentration_contours.rst

.. include:: post_processing/stress_contours.rst

.. include:: post_processing/first_invariant_of_strain.rst

.. include:: post_processing/second_invariant_of_strain.rst

.. include:: post_processing/third_invariant_of_strain.rst

.. include:: post_processing/velocity_divergence.rst

.. include:: post_processing/particle_velocity_divergence.rst

.. include:: post_processing/total_velocity_divergence.rst

.. include:: post_processing/electric_field.rst

.. include:: post_processing/electric_field_magnitude.rst

.. include:: post_processing/enormsq_field.rst

.. include:: post_processing/enormsq_field_norm.rst

.. include:: post_processing/viscosity.rst

.. include:: post_processing/density.rst

.. include:: post_processing/lame_mu.rst

.. include:: post_processing/lame_lambda.rst

.. include:: post_processing/von_mises_strain.rst

.. include:: post_processing/von_mises_stress.rst

.. include:: post_processing/navier_stokes_residuals.rst

.. include:: post_processing/moving_mesh_residuals.rst

.. include:: post_processing/mass_diffusion_vectors.rst

.. include:: post_processing/diffusive_mass_flux_vectors.rst

.. include:: post_processing/mass_fluxlines.rst

.. include:: post_processing/energy_conduction_vectors.rst

.. include:: post_processing/energy_fluxlines.rst

.. include:: post_processing/time_derivatives.rst

.. include:: post_processing/mesh_stress_tensor.rst

.. include:: post_processing/real_solid_stress_tensor.rst

.. include:: post_processing/mesh_strain_tensor.rst

.. include:: post_processing/viscoplastic_def_grad_tensor.rst

.. include:: post_processing/lagrangian_convection.rst

.. include:: post_processing/normal_and_tangent_vectors.rst

.. include:: post_processing/error_zz_velocity.rst

.. include:: post_processing/error_zz_heat_flux.rst

.. include:: post_processing/error_zz_pressure.rst

.. include:: post_processing/user_defined_post_processing.rst

.. include:: post_processing/porous_saturation.rst

.. include:: post_processing/total_density_of_solvents_in_porous_media.rst

.. include:: post_processing/density_of_solvents_in_gas_phase_in_porous_media.rst

.. include:: post_processing/density_of_liquid_phase_in_porous_media.rst

.. include:: post_processing/gas_phase_darcy_velocity_in_porous_media.rst

.. include:: post_processing/liquid_phase_darcy_velocity_in_porous_media.rst

.. include:: post_processing/capillary_pressure_in_porous_media.rst

.. include:: post_processing/grid_peclet_number_in_porous_media.rst

.. include:: post_processing/supg_velocity_in_porous_media.rst

.. include:: post_processing/vorticity_vector.rst
