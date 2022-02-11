Post Processing Fluxes and Data
###################################

By *Post-processing Fluxes* we mean area-integrated fluxes that can be calculated for any flux
quantity across any surface demarcated by a side set. The area-integrated flux is in fact a total
flow rate across the boundary. Examples include heat-flow, total force of a liquid on a surface, and species flow (both diffusive and convective). The integrated flux quantities are output to a
specified file at each time step, together with the time stamp and the convective and diffusive
components. This capability is useful for extracting engineering results from an analysis, and can
further be used to as an objective function evaluator for engineering optimization problems (cf.
*Post Processing Flux Sensitivities* card below).

*Post Processing Data* output can be used to produce spatial {*value, x, y, z*} sets on a specified side set of any primitive variable in the problem, viz. pressure, x-component of velocity, etc. The quantity *value* is the value of the variable at a node in the side set, and x, y, z are the coordinates of the node.

.. include:: post_processing_fluxes/post_processing_fluxes.rst

.. include:: post_processing_fluxes/flux.rst

.. include:: post_processing_fluxes/end_of_flux.rst

.. include:: post_processing_fluxes/post_processing_data.rst

.. include:: post_processing_fluxes/data.rst

.. include:: post_processing_fluxes/end_of_data.rst

.. include:: post_processing_fluxes/post_processing_flux_sensitivities.rst

.. include:: post_processing_fluxes/flux_sens.rst

.. include:: post_processing_fluxes/end_of_flux_sens.rst

.. include:: post_processing_fluxes/post_processing_date_sensitivities.rst

.. include:: post_processing_fluxes/data_sens.rst

.. include:: post_processing_fluxes/end_of_data_sens.rst
