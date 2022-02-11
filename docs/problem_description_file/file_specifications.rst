File Specification
######################
 
In general, this first section of the main input file is used to direct *Goma* I/O through a series of
named external files that contain information about the finite element mesh, the initial guess of a
solution vector, and output options for saving solutions for continuation, remesh, etc. The
required and optional input records are as follows:
 
.. include:: file_specifications/fem_file.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/output_exodusII_file.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/guess_file.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/soln_file.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/brk_file.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/write_intermediate_results.rst

-------------------------------------------------------------------------------

.. include:: file_specifications/write_initial_solution.rst
