3.2 Continuation Specifications
===============================

This section of input records is used to direct all automatic continuation procedures. The entire section is completely optional. Basically, automatic continuation can be accomplished in steady state simulations (see Time Integration card) through any one or combination of parameters. These parameters can be any one or combination of the input floats required on the boundary condition cards (see Section 4.10) or material property cards (see Chapter 5). The cards in this section are used to specify the parameters that will be marched automatically, the method of marching (e.g. zero-order, first-order, multiparameter first-order, etc.), the limits of parameter values, and other sundry options. Much of this capability can now be managed from the LOCA library package (Library of Continuation Algorithms - Salinger, et al. 2002).

.. toctree::
   :maxdepth: 2

   continuation
   continuation_type
   number_of_user_continuation_functions
   boundary_condition_id
   boundary_condition_data_float_tag
   material_id
   material_property_tag
   material_property_tag_subindex
   initial_parameter_value
   final_parameter_value
   delta_s
   maximum_number_of_path_steps
   maximum_path_value
   minimum_path_step
   maximum_path_step
   path_step_parameter
   path_step_error
   continuation_printing_frequency
   second_frequency
   loca_method
