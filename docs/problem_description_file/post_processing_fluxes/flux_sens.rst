*************
**FLUX_SENS**
*************

::

   FLUX_SENS = {flux_type} <bc_id> <blk_id> <species_id> {sens_type}
   <sens_id> <sens_flt> <file_name>

-----------------------
**Description / Usage**
-----------------------

FLUX_SENS cards request the calculation of the sensitivity of an integrated flux with
respect to a boundary condition or material parameter. As many FLUX_SENS cards as
desired can be input to *Goma*. Definitions of the input parameters are as follows:

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
|             |                                                                        |
|             |For every request, the specified sensitivity of the integrated          |
|             |diffusive flux followed by that of the convective portion               |
|             |over the requested boundary will be appended to the                     |
|             |specified file. If the convective flux is not applicable (cf.           |
|             |FLUX card), the second quantity will be zero. In all cases              |
|             |the area of the face (covered by the entire side set) and the           |
|             |time value are also output.                                             |
+-------------+------------------------------------------------------------------------+
|<bc_id>      |The boundary flag identifier, an integer associated with the            |
|             |boundary location (node set in EXODUS II) in the problem                |
|             |domain on which the integrated flux sensitivity is desired.             |
+-------------+------------------------------------------------------------------------+
|<blk_id>     |An integer that designates the mesh block (material) from               |
|             |which the flux sensitivity integral should be performed.                |
|             |This has implications on internal boundaries.                           |
+-------------+------------------------------------------------------------------------+
|<species_id> |An integer that identifies the species number if an integrated          |
|             |species flux sensitivity is requested.                                  |
+-------------+------------------------------------------------------------------------+
|{sens_type}  |A two-character entry that identifies the sensitivity type, where:      |
|             |                                                                        |
|             | * **BC** denotes a sensitivity w.r.t. to a boundary condition          |
|             |   parameter.                                                           |
|             | * **MT** denotes a sensitivity w.r.t. to a material property           |
|             |   parameter.                                                           |
+-------------+------------------------------------------------------------------------+
|<sens_id>    |An integer that identifies the sensitivity. If **BC** is specified      |
|             |for {sens_type}, then this value is the BC card number. If              |
|             |**MT** is specified for {sens_type}, this value is the material number. |
+-------------+------------------------------------------------------------------------+
|<sens_flt>   |A floating point number whose meaning is also dependent                 |
|             |on the selection of {sens_type}. If **BC** is specified, this           |
|             |value is the BC data float number. If **MT** is specified, this         |
|             |value is the material property tag.                                     |
+-------------+------------------------------------------------------------------------+
|<file_name>  |A character string corresponding to a file name into which              |
|             |these fluxes should be printed.                                         |
+-------------+------------------------------------------------------------------------+

------------
**Examples**
------------

Here is an example input deck:
::

   Post Processing Flux Sensitivities =
   FLUX_SENS = VOLUME_FLUX 1 1 0 BC 5 3 flux_sens.out
   END OF FLUX_SENS

-------------------------
**Technical Discussion**
-------------------------

Currently, the flux sensitivities do not account for the implicit sensitivity of material
properties. That is, :math:`dFORCE_X` / :math:`dBC_{float}` does not include a contribution from
:math:`d\mu` / :math:`dBC_{float}` , but :math:`dFORCE_X` / :math:`dMT_{property}`
should be done correctly. In addition,
sensitivities of integrated fluxes in solid materials have not been implemented yet.

NOTE: In order to compute flux sensitivities with respect to Dirichlet boundary
condition floats, the boundary conditions need to use the *residual method* in the input
file as described in the *Boundary Condition Specification* introduction, i.e. the optional
parameter should be set to 1.0.



--------------
**References**
--------------

No References.