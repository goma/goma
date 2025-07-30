Eigensolver Specifications
==========================

The ability to solve for the stability of a base flow is a very powerful tool. Often, the important 
characteristics of a flow can be summarized in the answer to the question "is the flow stable?". 
Although the following cards are in active use at the time of this writing, sweeping changes are 
coming to the eigensolver sections of Goma. In particular, the old code (called "eggroll") is being 
replaced with newer methods (in the ARPACK library), as well as being coupled to the 
continuation and tracking algorithms (in the LOCA library).

Input specifications for this section of input records is discussed in a separate, comprehensive 
manual (Gates, et. al., 2000); an update to this manual will be completed during the fall of 2002 
(Labreche, et. al., 2002). Either of these manuals contains a thorough discussion of how to 
successfully compute the stability and interesting modes of an underlying base flow.

.. toctree::
   :maxdepth: 1

   eigen_algorithm
   eigen_number_of_modes
   eigen_record_modes
   eigen_size_of_krylov_subspace
   eigen_maximum_iterations
   eigen_number_of_filter_steps
   eigen_recycle
   eigen_tolerance
   eigen_matrix_output
   eigen_initial_vector_weight
   eigen_initial_shifts
   eigen_wave_numbers
   eigen_cayley_sigma
   eigen_cayley_mu
   eigen_relative_tolerance
   eigen_linear_solver_tolerance
   eigenvalue_output_frequency
   eigenvector_output_frequency
   eigenvector_output_file
