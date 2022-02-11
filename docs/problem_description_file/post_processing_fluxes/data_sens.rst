*************
**DATA_SENS**
*************

::

   DATA = {data_type} <bc_id> <blk_id> <species_id> {sens_type} <sens_id>
   <sens_flt> <file_name>

-----------------------
**Description / Usage**
-----------------------

As many of these *DATA_SENS* cards as desired can be input to direct *Goma* to print the
sensitivity of a specified variable with respect to a boundary condition or material
parameter on a specified node set. Definitions of the input parameters are as follows:

+-------------+------------------------------------------------------------------------+
|{data_type}  |A keyword that can have any one of the following values:                |
|             |                                                                        |
|             | * **VELOCITY[1-3]**                                                    |
|             | * **TEMPERATURE**                                                      |
|             | * **MASS_FRACTION**                                                    |
|             | * **MESH_DISPLACEMENT[1-3]**                                           |
|             | * **SURFACE**                                                          |
|             | * **PRESSURE**                                                         |
|             | * **POLYMER_STRESS[1-3][1-3]**                                         |
|             | * **SOLID_DISPLACEMENT[1-3]**                                          |
|             | * **VELOCITY_GRADIENT[1-3][1-3]**                                      |
|             | * **VOLTAGE**                                                          |
|             | * **FILL**                                                             |
|             | * **SHEAR_RATE**                                                       |
|             | * **PVELOCITY[1-3]**                                                   |
|             | * **VolFracPh_[0-4]**                                                  |
|             | * **POR_LIQ_PRES**                                                     |
|             | * **POR_GAS_PRES**                                                     |
|             | * **POR_PORSITY**                                                      |
|             | * **POR_SATURATION**                                                   |
|             | * **VORT_DIR[1-3]**                                                    |
|             | * **VORT_LAMBDA**                                                      |
|             |                                                                        |
|             |Each request will result in the point coordinates and the               |
|             |specified sensitivity value being printed to the specified file.        |
+-------------+------------------------------------------------------------------------+
|<bc_id>      |The boundary flag identifier, an integer associated with the            |
|             |boundary location (node set in EXODUS II) in the problem                |
|             |domain on which the variable sensitivity is desired.                    |
+-------------+------------------------------------------------------------------------+
|<blk_id>     |An integer that designates the mesh block (material) from               |
|             |which the variable sensitivity should be taken. This has                |
|             |implications for discontinuous variables on internal boundaries.        |
+-------------+------------------------------------------------------------------------+
|<species_id> |An integer that identifies the species number if a species              |
|             |sensitivity is requested.                                               |
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
|             |the data should be printed.                                             |
+-------------+------------------------------------------------------------------------+

------------
**Examples**
------------

The following example shows a sample input deck section with three data requests:
::

   Post Processing Data Sensitivities =
   DATA_SENS = VELOCITY2 1 1 0 BC 5 3 data_sens.out
   DATA_SENS = VELOCITY1 6 1 0 BC 5 3 data_sens.out
   DATA_SENS = VELOCITY1 6 1 0 BC 4 0 data_sens.out
   END OF DATA_SENS

-------------------------
**Technical Discussion**
-------------------------

NOTE: In order to compute data sensitivities with respect to dirichlet boundary
condition floats, the boundary conditions need to be "soft" set in the input file, i.e. the
optional parameter should be set to 1.0.



--------------
**References**
--------------

No References.