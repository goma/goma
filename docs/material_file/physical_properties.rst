Physical Properties
#######################

The intrinsic property of materials essential to *Goma* is the density. As *Goma* presumes that all
materials are incompressible, density is a constant in the governing differential equations.
However, several options for models of density are present in the code because numerous
processes lead to density changes, though during any analysis cycle, the density is constant.

.. toctree::
   :maxdepth: 1

   physical_properties/default_database
   physical_properties/density
