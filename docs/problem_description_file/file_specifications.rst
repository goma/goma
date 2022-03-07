File Specification
######################
 
In general, this first section of the main input file is used to direct *Goma* I/O through a series of
named external files that contain information about the finite element mesh, the initial guess of a
solution vector, and output options for saving solutions for continuation, remesh, etc. The
required and optional input records are as follows:
 
.. toctree::
   :maxdepth: 1

   file_specifications/fem_file
   file_specifications/output_exodusII_file
   file_specifications/guess_file
   file_specifications/soln_file
   file_specifications/write_intermediate_results
   file_specifications/write_initial_solution
   file_specifications/external_decomposition
   file_specifications/decomposition_type
