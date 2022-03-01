Time Integration Specifications
###################################

The first card in this section dictates whether the problem is a steady state or transient simulation.
This card is required. If the *steady state* option is chosen, then the remaining input records are not
required, as the rest of the records are used to set parameters for transient simulations, e.g., time
step size, time step error control, etc. Some records are optional even for a transient simulation, as
indicated below. It should be noted that the mass-matrix term multiplier in the *Problem
Description section* (see, for example, the *EQ=* cards), must be set to one (1) for the transient run
to evolve the fields in time. The only equations that are taken as purely quasi static are the
*EQ=mesh* equations for the situation in which the *Mesh Motion* type is *Arbitrary*.

In addition to the transient parameter information, some Level-Set function information is also
supplied to *Goma* in this section. The method of Level-Sets is used to track fluid-fluid or fluidsolid
interfaces in an Eulerian fashion, making the problem inherently transient.

.. toctree::
   :maxdepth: 1

   time_integration/time_integration
   time_integration/delta_t
   time_integration/maximum_number_of_time_steps
   time_integration/maximum_time
   time_integration/minimum_time_step
   time_integration/maximum_time_step
   time_integration/minimum_resolved_time_step
   time_integration/courant_number_limit
   time_integration/time_step_parameter
   time_integration/time_step_error
   time_integration/printing_frequency
   time_integration/fix_frequency
   time_integration/second_frequency_time
   time_integration/initial_time

