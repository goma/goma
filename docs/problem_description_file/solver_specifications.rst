Solver Specifications
#########################

This required section directs the nonlinear iteration strategy with associated parameters (e.g.,
Newton’s method options), matrix solution strategy and parameters, and other sundry options and
toggles for the pressure stabilization approach and linear stability analysis capability. With regard
to the parameters associated with matrix solution methods, it is important to understand that there
are two major classes of solvers - direct and iterative solvers. Direct solvers are the most robust,
but can be computationally impractical for some larger systems. Iterative solvers and associated
preconditioners are the only practical options for large-scale problems (viz., very large twodimensional
problems and virtually all three-dimensional problems). Choosing the solver settings
for good convergence of iterative matrix solvers can be an artful task for Navier-Stokes problems
and other poorly conditioned systems. It is recommended that the user consult the comprehensive
report by Schunk, et al. (2002) for an overview and further usage tips.

.. include:: solver_specifications/solution_algorithm.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_storage_format.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/stratimikos_file.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/preconditioner.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_subdomain_solver.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_scaling.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_residual_norm_type.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_output_type.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_factorization_reuse.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_graph_fillin.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_factorization_overlap.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_overlap_type.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_auxiliary_vector.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_drop_tolerance.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_polynomial_order.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_reorder.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_factorization_save.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_ilut_fill_factor.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_rilu_relax_factor.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_bilu_threshold.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_relative_threshold.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/matrix_absolute_threshold.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/size_of_krylov_subspace.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/orthogonalization.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/maximum_linear_solve_iterations.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/number_of_newton_iterations.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/modified_newton_tolerance.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/jacobian_reform_time_stride.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/newton_correction_factor.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/normalized_residual_tolerance.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/normalized_correction_tolerance.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/residual_ratio_tolerance.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/pressure_stabilization.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/pressure_stabilization_scaling.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/linear_stability.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/filter_concentration.rst

-------------------------------------------------------------------------------

.. include:: solver_specifications/disable_viscosity_sensitivities.rst

