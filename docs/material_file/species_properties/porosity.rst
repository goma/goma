************
**Porosity**
************

::

   Porosity = {model_name} <float1> [float2]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the porosity model for the anode or separator or cathode
region in a thermal battery cell.

Definitions of the {model_name} and the associated input parameters (<float>) are as
follows:

+----------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**          |the name of the porosity model.                                                      |
|                      |                                                                                     |
|                      | * {float1} - the porosity value.                                                    |
+----------------------+-------------------------------------------------------------------------------------+
|**THERMAL_BATTERY**   |the name of the porosity model.                                                      |
|                      |                                                                                     |
|                      | * <float1> - the initial value of porosity                                          |
|                      | * <float2> - specifies the change of molar volume in the anode or cathode electrode |
|                      |   material per electron transferred, as stated in                                   |
|                      |                                                                                     |
|                      |   .. figure:: /figures/430_goma_physics.png                                         |
|                      |      :align: center                                                                 |
|                      |      :width: 90%                                                                    |
|                      |                                                                                     |
|                      |   where si is stoichiometric coefficient of species or phase i, V is molar volume of|
|                      |   species or phase i, n is the number of electrons transfer in the anodic or        |
|                      |   cathodic electrochemical reaction, and the summation is over the number of solid  |
|                      |   phases.                                                                           |
+----------------------+-------------------------------------------------------------------------------------+

------------
**Examples**
------------

A sample input card for this material property might look like this:

::

   Porosity = THERMAL_BATTERY 0.244 8.1185

-------------------------
**Technical Discussion**
-------------------------

* This is a porosity model for a special application in which the model for the
  diffusion constitutive equation is *STEFAN_MAXWELL_CHARGED*, which
  enables modeling the transport of multiple charged species with simultaneous
  electrochemical reaction(s) in a concentrated solution, as in a thermal-battery cell.

* See the reference below for a discussion of Thermal Battery modeling with *Goma*.



--------------
**References**
--------------

SAND2000-0207: Final Report on LDRD Project: A Phenomenological Model for
Multicomponent Transport with Simultaneous Electrochemical Reactions in
Concentrated Solutions, K. S. Chen, G. H. Evans, R. S. Larson, D. R. Noble and W. G.
Houf, January 2000.