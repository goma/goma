
Problem Description
#######################

This section directs all input specifications required for differential equations, material type, mesh motion type, coordinate system, finite element basis function type, and several other input tasks. This section of input records, with the exception of the *Number of Materials* card (the first one
below), must be repeated for each material region in the problem. Within that region of the problem domain (and the corresponding section of the input file) there are no restrictions as to which differential or constraint equations can be specified, which is a unique capability of *Goma*.

However, some combinations or specifications do not make much sense, e.g., a cylindrical coordinate region combined with a cartesian one. It is recommended that the user consult the usage tutorials and example problems to get a feel for how this section is constructed.

.. toctree::
   :maxdepth: 1

   problem_description/number_of_materials
   problem_description/mat
   problem_description/coordinate_system
   problem_description/element_mapping
   problem_description/mesh_motion
   problem_description/number_of_bulk_species
   problem_description/material_is_nondilute
   problem_description/number_of_bulk_species_equations
   problem_description/default_material_species_type
   problem_description/number_of_viscoelastic_modes
   problem_description/number_of_matrices
   problem_description/matrix
   problem_description/disable_time_step_control
   problem_description/normalized_residual_tolerance
   problem_description/number_of_eq
   problem_description/energy
   problem_description/momentum
   problem_description/pmomentum
   problem_description/stress
   problem_description/eddy_visc
   problem_description/species_bulk
   problem_description/mesh
   problem_description/mom_solid
   problem_description/continuity
   problem_description/fill
   problem_description/lagr_mult_1_lagr_mult_2_lagr_mult_3
   problem_description/level_set
   problem_description/voltage
   problem_description/efield
   problem_description/enorm
   problem_description/shear_rate
   problem_description/vort_dir
   problem_description/vort_lambda
   problem_description/porous_sat
   problem_description/porous_unsat
   problem_description/porous_liq
   problem_description/porous_gas
   problem_description/porous_deform
   problem_description/porous_energy
   problem_description/surf_charge
   problem_description/shell_tension
   problem_description/shell_curvature
   problem_description/shell_angle
   problem_description/shell_diff_flux
   problem_description/shell_diff_curv
   problem_description/shell_normal
   problem_description/shell_surf_curv
   problem_description/shell_surf_div_v
   problem_description/grad_v_dot_n1_grad_v_dot_n2_grad_v_dot_n3
   problem_description/n_dot_curl_v
   problem_description/acous_preal
   problem_description/acous_pimag
   problem_description/acous_reyn_stress
   problem_description/potential1
   problem_description/potential2
   problem_description/lubp
   problem_description/lubp_2
   problem_description/shell_energy
   problem_description/shell_filmp
   problem_description/shell_filmh
   problem_description/shell_partc
   problem_description/shell_sat_closed
   problem_description/shell_sat_gasn
   problem_description/shell_sat_open
   problem_description/shell_sat_open_2
   problem_description/shell_deltah
   problem_description/moment
   problem_description/end_of_eq
   problem_description/end_of_mat
