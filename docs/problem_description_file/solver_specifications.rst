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

.. toctree::
   :maxdepth: 1

   solver_specifications/total_number_of_matrices
   solver_specifications/solution_algorithm
   solver_specifications/matrix_storage_format
   solver_specifications/stratimikos_file
   solver_specifications/preconditioner
   solver_specifications/matrix_subdomain_solver
   solver_specifications/matrix_scaling
   solver_specifications/matrix_residual_norm_type
   solver_specifications/matrix_output_type
   solver_specifications/matrix_factorization_reuse
   solver_specifications/matrix_graph_fillin
   solver_specifications/matrix_factorization_overlap
   solver_specifications/matrix_overlap_type
   solver_specifications/matrix_auxiliary_vector
   solver_specifications/matrix_drop_tolerance
   solver_specifications/matrix_polynomial_order
   solver_specifications/matrix_reorder
   solver_specifications/matrix_factorization_save
   solver_specifications/matrix_ilut_fill_factor
   solver_specifications/matrix_rilu_relax_factor
   solver_specifications/matrix_bilu_threshold
   solver_specifications/matrix_relative_threshold
   solver_specifications/matrix_absolute_threshold
   solver_specifications/size_of_krylov_subspace
   solver_specifications/orthogonalization
   solver_specifications/maximum_linear_solve_iterations
   solver_specifications/number_of_newton_iterations
   solver_specifications/modified_newton_tolerance
   solver_specifications/jacobian_reform_time_stride
   solver_specifications/newton_correction_factor
   solver_specifications/normalized_residual_tolerance
   solver_specifications/normalized_correction_tolerance
   solver_specifications/residual_ratio_tolerance
   solver_specifications/pressure_stabilization
   solver_specifications/pressure_stabilization_scaling
   solver_specifications/linear_stability
   solver_specifications/filter_concentration
   solver_specifications/disable_viscosity_sensitivities

