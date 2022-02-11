
Problem Description
#######################

This section directs all input specifications required for differential equations, material type, mesh motion type, coordinate system, finite element basis function type, and several other input tasks. This section of input records, with the exception of the *Number of Materials* card (the first one
below), must be repeated for each material region in the problem. Within that region of the problem domain (and the corresponding section of the input file) there are no restrictions as to which differential or constraint equations can be specified, which is a unique capability of *Goma*.

However, some combinations or specifications do not make much sense, e.g., a cylindrical coordinate region combined with a cartesian one. It is recommended that the user consult the usage tutorials and example problems to get a feel for how this section is constructed.

.. include:: problem_description/number_of_materials.rst

.. include:: problem_description/mat.rst

.. include:: problem_description/coordinate_system.rst

.. include:: problem_description/element_mapping.rst

.. include:: problem_description/mesh_motion.rst

.. include:: problem_description/number_of_bulk_species.rst

.. include:: problem_description/material_is_nondilute.rst

.. include:: problem_description/number_of_bulk_species_equations.rst

.. include:: problem_description/default_material_species_type.rst

.. include:: problem_description/number_of_viscoelastic_modes.rst

.. include:: problem_description/number_of_eq.rst

.. include:: problem_description/energy.rst

.. include:: problem_description/momentum.rst

.. include:: problem_description/pmomentum.rst

.. include:: problem_description/stress.rst

.. include:: problem_description/species_bulk.rst

.. include:: problem_description/mesh.rst

.. include:: problem_description/mom_solid.rst

.. include:: problem_description/continuity.rst

.. include:: problem_description/fill.rst

.. include:: problem_description/lagr_mult_1_lagr_mult_2_lagr_mult_3.rst

.. include:: problem_description/level_set.rst

.. include:: problem_description/voltage.rst

.. include:: problem_description/efield.rst

.. include:: problem_description/enorm.rst

.. include:: problem_description/shear_rate.rst

.. include:: problem_description/vort_dir.rst

.. include:: problem_description/vort_lambda.rst

.. include:: problem_description/porous_sat.rst

.. include:: problem_description/porous_unsat.rst

.. include:: problem_description/porous_liq.rst

.. include:: problem_description/porous_gas.rst

.. include:: problem_description/porous_deform.rst

.. include:: problem_description/porous_energy.rst

.. include:: problem_description/surf_charge.rst

.. include:: problem_description/shell_tension.rst

.. include:: problem_description/shell_curvature.rst

.. include:: problem_description/shell_angle.rst

.. include:: problem_description/shell_diff_flux.rst

.. include:: problem_description/shell_diff_curv.rst

.. include:: problem_description/shell_normal.rst

.. include:: problem_description/shell_surf_curv.rst

.. include:: problem_description/shell_surf_div_v.rst

.. include:: problem_description/grad_v_dot_n1_grad_v_dot_n2_grad_v_dot_n3.rst

.. include:: problem_description/n_dot_curl_v.rst

.. include:: problem_description/acous_preal.rst

.. include:: problem_description/acous_pimag.rst

.. include:: problem_description/acous_reyn_stress.rst

.. include:: problem_description/potential1.rst

.. include:: problem_description/potential2.rst

.. include:: problem_description/lubp.rst

.. include:: problem_description/lubp_2.rst

.. include:: problem_description/shell_energy.rst

.. include:: problem_description/shell_filmp.rst

.. include:: problem_description/shell_filmh.rst

.. include:: problem_description/shell_partc.rst

.. include:: problem_description/shell_sat_closed.rst

.. include:: problem_description/shell_sat_gasn.rst

.. include:: problem_description/shell_sat_open.rst

.. include:: problem_description/shell_sat_open_2.rst

.. include:: problem_description/shell_deltah.rst

.. include:: problem_description/end_of_eq.rst

.. include:: problem_description/end_of_mat.rst
