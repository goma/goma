**************
**VOLUME_INT**
**************

::

   VOLUME_INT = {volume_type} <blk_id> <species_no> <file_string> [float_list]

-----------------------
**Description / Usage**
-----------------------

The VOLUME_INT card activates computation of specified volumetric integrals
during post processing. As many of these VOLUME_INT cards as desired can be input
to *Goma*. Definitions of the input parameters are as follows:

+----------------+----------------------------------------------------------------+
|<volume_type>   |Several choices of volumetric integral are allowed and are      |
|                |referenced through this parameter. The permissible values       |
|                |and corresponding volume integral follow:                       |
|                |                                                                |
|                | * **VOLUME**-Volume of the element block specified by <blk_id>.|
|                | * **DISSIPATION**-Total viscous dissipation, :math:`\tau`:     |
|                |   :math:`\Delta v` , in the element block specified            |
|                |   by <blk_id>.                                                 |
|                | * **JOULE**-Total Joule or Ohmic heating,                      |
|                |   :math:`\frac{1}{\sigma}` (:math:`\underline{J}` •            |
|                |   :math:`\underline{J}`) , in the element block specified      |
|                |   by <blk_id>.                                                 |
|                | * **SPECIES_MASS**-Integral of concentration of the            |
|                |   component specified by <species_no> in the element           |
|                |   block specified by < blk_id>.                                |
|                | * **MOMENTUMX, MOMENTUMY,** or **MOMENTUMZ**-Integral of       |
|                |   appropriate component of the momentum flux                   |
|                |   :math:`\rho\vec{v}` over the element block <blk_id>.         |
|                | * **STRESS_TRACE**-Integral of the trace of the complete       |
|                |   stress tensor (-:math:`p\zeta + \tau`) over the element block|
|                |   *blk_id*.                                                    |
|                | * **HEAT_ENERGY**-Integral of the sensible heat over           |
|                |   <blk_id> ( not currently implemented).                       |
|                | * **POSITIVE_FILL, NEGATIVE_FILL**-Volume                      |
|                |   integral of region occupied by positive (negative)           |
|                |   values of the FILL variable in element block < blk_id>.      |
|                |   Note, for either of these cards, [float_list] is required.   |
|                |   **NOTE for Level-Set users**: There are numerous other       |
|                |   quantities (too-lengthy and esoteric to list here) that can  |
|                |   be integraded vis-a-vis level set fields. Please see code.   |
|                | * **NEGATIVE_VX, NEGATIVE_VY, NEGATIVE_VZ**-Velocity integral  |
|                |   in one of the three directions over just                     |
|                |   the region occupied by negative values of the FILL           |
|                |   variable in level set problems. Note, for any of these       |
|                |   cards, the [float_list] is required.                         |
|                | * **POROUS_LIQ_INVENTORY**-Volume integral of                  |
|                |   bulk liquid component density (gas and liquid phase) in      |
|                |   a porous medium. Result is a total inventory of liquid       |
|                |   in the porous phase.                                         |
|                | * **SPEED_SQUARED**-Volume integral of the square of           |
|                |   the speed, viz. :math:`\underline{v}` • :math:`\underline{v}`|
|                |   . Used                                                       |
|                |   to measure norm of fluid kinetic energy level. Should tend   |
|                |   to zero for a fluid at rest.                                 |
|                | * **USER**-Volume integral is supplied by the user (not        |
|                |   currently implemented).                                      |
|                | * **SURFACE_SPECIES**-Generate locus of points which           |
|                |   correspond to a surface of constant species                  |
|                |   concentration according to :math:`Ac_1` + :math:`Bc_2` +     |
|                |   :math:`Cc_3` + D=0.                                          |
|                |   Currently only implemented for 3D linear elements.           |
|                | * **LUB_LOAD**-“Volume integral” of lubrication pressure       |
|                |   over entire mesh shell block, which is useful for            |
|                |   computing the overall lubrication load. This is              |
|                |   actually an area integral over the shell, thereby yielding   |
|                |   a force.                                                     |
|                | * **ELOADX; ELOADY; ELOADZ**-Volume integral of                |
|                |   electric field or the gradient of the electric potential for |
|                |   electrostatic problems.                                      |
|                | * **RATE_OF_DEF_II**-Volume integral of the second             |
|                |   invariant of the rate-of-deformation tensor.                 | 
+----------------+----------------------------------------------------------------+
|<blk_id>        |The element block id for which the volume integral is requested.|
+----------------+----------------------------------------------------------------+
|<species_no>    |The species number for **SPECIES_MASS** volume integral.        |
+----------------+----------------------------------------------------------------+
|<file_string>   |A character string that corresponds to the name of the text     |
|                |file that will receive the results of the integration at each   |
|                |time step.                                                      |
+----------------+----------------------------------------------------------------+
|[float_list]    |A floating point value that specifies the length scale of the   |
|                |smooth Heaviside function. This parameter is only used for      |
|                |*VOLUME_INT* cards in which the {volume_type} is                |
|                |**{POSITIVE|NEGATIVE} _FILL** or                                |
|                |**NEGATIVE_V{X|Y|Z}**. The float list is also used for the      |
|                |constants A, B, C, etc in the **SURFACE_SPECIES** type.         |
+----------------+----------------------------------------------------------------+

------------
**Examples**
------------

Here is an example of an input deck with 3 VOLUME_INT cards.
::

   Post Processing Volumetric Integration =
   VOLUME_INT = VOLUME 1 0 volume.out
   VOLUME_INT = SPECIES_MASS 2 3 species3.out
   VOLUME_INT = NEGATIVE_FILL 1 0 fill.out 0.1
   END OF VOLUME_INT

-------------------------
**Technical Discussion**
-------------------------

The volume integrations are carried out as follows:

+------------------------+---------------------------------------------------------------+
|**volume_type**         |**volume integral**                                            |
+------------------------+---------------------------------------------------------------+
|VOLUME                  |:math:`\int` dV                                                |
+------------------------+---------------------------------------------------------------+
|DISSIPATION             |:math:`\int` (-p :math:`\zeta + \tau`) • :math:`\Delta` vdV    |
+------------------------+---------------------------------------------------------------+
|JOULE                   |:math:`\int` :math:`\frac{1}{\sigma}` J • JdV                  |
+------------------------+---------------------------------------------------------------+
|SPECIES_MASS            |:math:`\int` :math:`c_jdV`                                     |
+------------------------+---------------------------------------------------------------+
|MOMENTUM_{X|Y|Z}        |:math:`\int` :math:`\rho` (i|j|k) • vdV                        |
+------------------------+---------------------------------------------------------------+
|STRESS_TRACE            |:math:`\int` tr(-p :math:`\zeta + \tau`) dV                    |
+------------------------+---------------------------------------------------------------+
|{POSITIVE|NEGATIVE}_FILL|:math:`\int` H(:math:`\phi`) dV                                |
+------------------------+---------------------------------------------------------------+
|NEGATIVE_V{X|Y|Z}       |:math:`\int` H(:math:`\phi`) {i|j|k} • vdV                     |
+------------------------+---------------------------------------------------------------+
|POROUS_LIQ_INVENTORY    |:math:`\int` [:math:`\rho_{gas}` :math:`\phi` (1-S) +          |
|                        |:math:`\rho_{liq}` :math:`\phi` S] dV                          |
+------------------------+---------------------------------------------------------------+



--------------
**References**
--------------

No References.