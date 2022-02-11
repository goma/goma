********
**FLUX**
********

::

   FLUX = {flux_type} <bc_id> <blk_id> <species_id> <file_name> [profile]

-----------------------
**Description / Usage**
-----------------------

FLUX cards are used to calculate the integrated fluxes of momentum, mass, energy,
etc. on a specified side set during post processing. As many of these FLUX cards as
desired can be input to *Goma* to direct the calculations. For example, multiple cards
may be used to compute a particular flux, e.g. FORCE_NORMAL, on different side
sets or different fluxes on the same side set. Cards with identical fluxes and identical
side sets could be used to output the flux calculations to different files. Definitions of
the input parameters are:

+-------------+------------------------------------------------------------------------+
|{flux_type}  |A keyword that can have any one of the following values:                |
|             |                                                                        |
|             | * **FORCE_NORMAL**                                                     |
|             | * **FORCE_TANGENT1**                                                   |
|             | * **FORCE_TANGENT2**                                                   |
|             | * **FORCE_X**                                                          |
|             | * **FORCE_Y**                                                          |
|             | * **FORCE_Z**                                                          |
|             | * **VOLUME_FLUX**                                                      |
|             | * **SPECIES_FLUX**                                                     |
|             | * **HEAT_FLUX**                                                        |
|             | * **TORQUE**                                                           |
|             | * **AVERAGE_CONC**                                                     |
|             | * **SURF_DISSIP**                                                      |
|             | * **AREA**                                                             |
|             | * **VOL_REVOLUTION**                                                   |
|             | * **PORE_LIQ_FLUX**                                                    |
|             | * **CHARGED_SPECIES_FLUX**                                             |
|             | * **CURRENT_FICKIAN**                                                  |
|             | * **CURRENT**                                                          |
|             | * **ELEC_FORCE_NORMAL**                                                |
|             | * **ELEC_FORCE_TANGENT1**                                              |
|             | * **ELEC_FORCE_TANGENT2**                                              |
|             | * **ELEC_FORCE_X**                                                     |
|             | * **ELEC_FORCE_Y**                                                     |
|             | * **ELEC_FORCE_Z**                                                     |
|             | * **NET_SURF_CHARGE**                                                  |
|             | * **ACOUSTIC_FLUX_NORMAL**                                             |
|             | * **ACOUSTIC_FLUX_TANGENT1**                                           |
|             | * **ACOUSTIC_FLUX_TANGENT2**                                           |
|             | * **ACOUSTIC_FLUX_X**                                                  |
|             | * **ACOUSTIC_FLUX_Y**                                                  |
|             | * **ACOUSTIC_FLUX_Z**                                                  |
|             |                                                                        |
|             |For every request, the integral of the diffusive portion                |
|             |followed by that of the convective portion over the                     |
|             |requested boundary will be appended to the specified file. If           |
|             |the convective flux is not applicable (i.e.for flux_types               |
|             |**VOLUME_FLUX, TORQUE, AVERAGE_CONC** and                               |
|             |**AREA**), the second quantity will be zero. In all cases the           |
|             |area of the face (covered by the entire side set) and the time          |
|             |value are also output.                                                  |
+-------------+------------------------------------------------------------------------+
|<bc_id>      |The boundary flag identifier, an integer associated with the            |
|             |boundary location (side set in EXODUS II) in the problem                |
|             |domain on which the integrated flux is desired.                         |
+-------------+------------------------------------------------------------------------+
|<blk_id>     |An integer that designates the mesh block (material) from               |
|             |which the flux integral should be performed. This has                   |
|             |implications on internal boundaries.                                    |
+-------------+------------------------------------------------------------------------+
|<species_id> |An integer that identifies the species number if an integrated          |
|             |species flux is requested.                                              |
+-------------+------------------------------------------------------------------------+
|<file_name>  |A character string corresponding to a file name into which              |
|             |these fluxes should be printed.                                         |
+-------------+------------------------------------------------------------------------+
|[profile]    |Inclusion of the optional string “profile’ to this card will            |
|             |cause the coordinates (x,y,z), the diffusive integrand, and             |
|             |the convective integrand at each integration point to be                |
|             |printed to the file designated above. You can, for example,             |
|             |print out a pressure distribution used to compute a force.              |
+-------------+------------------------------------------------------------------------+

------------
**Examples**
------------

The following example shows a sample input deck section that requests five such
integrated fluxes:
::

   Post Processing Fluxes =

::

   FLUX = FORCE_X 5 1 0 side5.out

::

   FLUX = FORCE_Y 5 1 0 side5prof.out   profile

::

   FLUX = FORCE_NORMAL 8 1 0 side8.out

::

   FLUX = FORCE_TANGENT1 8 1 0 side8.out

::

   FLUX = VOLUME_FLUX 8 1 0 side8.out

::

   END OF FLUX

-------------------------
**Technical Discussion**
-------------------------

The permissible flux types are those listed in file mm_post_def.h for struct
*Post_Processing_Flux_Names, pp_flux_names* being one variable of this struct type.

The flux integrations are carried out as follows:

+----------------------+--------------------------------------------+-----------------------------+
|**FLUX**              |**DIFFUSIVE FLUX**                          |**CONVECTIVE FLUX**          |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_NORMAL          |:math:`\int` n • :math:`\underline{T}` • ndA|:math:`\int\rho` n •         |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_TANGENT1        |:math:`\int_1` • :math:`\underline{T}` • ndA|:math:`\int\rho` :math:`t_1` |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_TANGENT2        |:math:`\int_2` • :math:`\underline{T}` • ndA|:math:`\int\rho` :math:`t_2` |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_X               |:math:`\int` i • :math:`\underline{T}` • ndA|:math:`\int\rho` i           |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_Y               |:math:`\int` j • :math:`\underline{T}` • ndA|:math:`\int\rho` j           |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|FORCE_Z               |:math:`\int` k • :math:`\underline{T}` • ndA|:math:`\int\rho` k           |
|                      |                                            |(v - :math:`v_m`) v • ndA    |
+----------------------+--------------------------------------------+-----------------------------+
|VOLUME_FLUX           |:math:`\int` n • (v - :math:`v_m`) dA       |for ARBITRARY mesh motion.   |
+----------------------+--------------------------------------------+-----------------------------+
|                      |:math:`\int` n • ddA                        |for LAGRANGIAN mesh motion.  |
+----------------------+--------------------------------------------+-----------------------------+
|SPECIES_FLUX          |:math:`\int` (-:math:`D_jn` •               |:math:`\int\rho` n •         |
|                      |:math:`\Delta` cj) dA                       |( v - :math:`v_m` ) cjdA     |
+----------------------+--------------------------------------------+-----------------------------+
|HEAT_FLUX             |:math:`\int` (-kn • :math:`\Delta` T) dA    |:math:`\int\rho` CpTn •      |
|                      |                                            |( v - :math:`v_m` ) dA       |
+----------------------+--------------------------------------------+-----------------------------+
|TORQUE                |:math:`\int` :math:`re_r` ×                 |                             |
|                      |( :math:`\underline{T}` • n) dA             |                             |
+----------------------+--------------------------------------------+-----------------------------+
|AVERAGE_CONC          |:math:`\int` cjdA                           |                             |
+----------------------+--------------------------------------------+-----------------------------+
|SURF_DISSIP           |:math:`\int\sigma\Delta` v •                |                             |
|                      |( :math:`\zeta` - nn) dA                    |                             |
+----------------------+--------------------------------------------+-----------------------------+
|AREA                  |:math:`\int` dA                             |                             |
+----------------------+--------------------------------------------+-----------------------------+
|VOL_REVOLUTION        |:math:`\int\frac{1}{2}`                     |                             |
|                      |:math:`\frac{r}{\sqrt{}{1 + (dr/dz)^2}}` dA |                             |
+----------------------+--------------------------------------------+-----------------------------+
|POR_LIQ_FLUX          |:math:`\int` n •                            |                             |
|                      |(:math:`\rho_lv_{darcy}`) dA                |                             |
+----------------------+--------------------------------------------+-----------------------------+
|CHARGED_SPECIES_FLUX  |:math:`\int` (-Djn • :math:`\Delta` cj) dA  |:math:`\int\rho` n •         |
|                      |                                            |( v - :math:`v_m` ) cjdA     |
+----------------------+--------------------------------------------+-----------------------------+
|CURRENT_FICKIAN       |:math:`\int` (-Djn • :math:`\Delta` cj) dA  |:math:`\int\rho` n •         |
|                      |                                            |( v - :math:`v_m` ) cjdA     |
+----------------------+--------------------------------------------+-----------------------------+
|PVELOCITY[1-3]        |:math:`\int` n • pvjdA                      |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_NORMAL     |:math:`\int` n :math:`\underline{T}_e` • ndA|                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_TANGENT1   |:math:`\int` :math:`t_1` •                  |                             |
|                      |:math:`\underline{T}_e` • ndA               |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_TANGENT2   |:math:`\int` :math:`t_2` •                  |                             |
|                      |:math:`\underline{T}_e` • ndA               |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_X          |:math:`\int` i •                            |                             |
|                      |:math:`\underline{T}_e` • ndA               |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_Y          |:math:`\int` j •                            |                             |
|                      |:math:`\underline{T}_e` • ndA               |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ELEC_FORCE_Y          |:math:`\int` k •                            |                             |
|                      |:math:`\underline{T}_e` • ndA               |                             |
+----------------------+--------------------------------------------+-----------------------------+
|NET_SURF_CHARGE       |:math:`\int` (-:math:`\varepsilon`          |                             |
|                      |:math:`\underline{n}` •                     |                             |
|                      |:math:`\underline{E}`) dA                   |                             |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_NORMAL  |:math:`\int` (-:math:`\frac{1}{kR}` n •     |                             |
|                      |:math:`\Delta P_{imag}`) dA                 |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}` n •   |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_TANGENT1|:math:`\int` (-:math:`\frac{1}{kR}`         |                             |
|                      |:math:`t_1` • :math:`\Delta P_{imag}`) dA   |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}`       |
|                      |                                            |:math:`t_1` •                |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_TANGENT2|:math:`\int` (-:math:`\frac{1}{kR}`         |                             |
|                      |:math:`t_2` • :math:`\Delta P_{imag}`) dA   |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}`       |
|                      |                                            |:math:`t_2` •                |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_X       |:math:`\int` (-:math:`\frac{1}{kR}`         |                             |
|                      |:math:`i` • :math:`\Delta P_{imag}`) dA     |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}`       |
|                      |                                            |:math:`i` •                  |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_Y       |:math:`\int` (-:math:`\frac{1}{kR}`         |                             |
|                      |:math:`j` • :math:`\Delta P_{imag}`) dA     |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}`       |
|                      |                                            |:math:`j` •                  |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+
|ACOUSTIC_FLUX_Z       |:math:`\int` (-:math:`\frac{1}{kR}`         |                             |
|                      |:math:`k` • :math:`\Delta P_{imag}`) dA     |                             |
|                      |                                            |:math:`\int`                 |
|                      |                                            |(-:math:`\frac{1}{kR}`       |
|                      |                                            |:math:`k` •                  |
|                      |                                            |:math:`\Delta P_{real}`) dA  |
+----------------------+--------------------------------------------+-----------------------------+

The SURF_DISSIP card is used to compute the energy dissipated at a surface by
surface tension (Batchelor, 1970). The VOL_REVOLUTION card is used in axi-
symmetric problems to compute the volume swept by revolving a surface around the
axis of symmetry (z-axis). Even though every flux card results in the area computation
of the side set, the AREA card is used when the area of a surface is part of an
augmenting condition. The POR_LIQ_FLUX term is valid only for saturated media
and the Darcy velocity is defined by :math:`\nu_{darcy}` = (:math:`\kappa` / :math:`\mu`) :math:`\Delta` 
:math:`p_{liq}` . For the more general case, refer to the *POROUS_LIQ_FLUX_CONST* boundary condition.



--------------
**References**
--------------

Batchelor, JFM, 1970. ..... need to fill-in reference; get from RBS

For information on using flux calculations as part of augmenting conditions, see:

SAND2000-2465: Advanced Capabilities in Goma 3.0 - Augmenting Conditions,
Automatic Continuation, and Linear Stability Analysis, I. D. Gates, I. D.,
Labreche, D. A. and Hopkins, M. M. (January 2001).