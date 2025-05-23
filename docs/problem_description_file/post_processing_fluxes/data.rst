********
**DATA**
********

::

   DATA = {data_type} <bc_id> <blk_id> <species_id> <file_name>

-----------------------
**Description / Usage**
-----------------------

A DATA card directs *Goma* to output the indicated primitive variable on a specified
node set. As many of these DATA cards as desired can be input to *Goma*. For example,
multiple cards may be used to output a particular variable, e.g. VELOCITY1, on
different node sets or different variables on the same node set. Cards with identical
variables and identical node sets could be used to output the variables to different files.
Definitions of the input parameters are as follows:

+-------------+------------------------------------------------------------------------+
|{data_type}  |A keyword that can have any one of the following primitive values:      |
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
|             | * **POLYMER_STRESS[1-3][1-3]_[1-7}**                                   |
|             | * **SPECIES_UNK[0-29]**                                                |
|             | * **VolFracPh_[0-4]**                                                  |
|             | * **POR_LIQ_PRES**                                                     |
|             | * **POR_GAS_PRES**                                                     |
|             | * **POR_PORSITY**                                                      |
|             | * **POR_SATURATION**                                                   |
|             | * **VORT_DIR[1-3]**                                                    |
|             | * **VORT_LAMBDA**                                                      |
|             |                                                                        |
|             |Each request will result in the point coordinates and the               |
|             |quantity value being printed to the specified file.                     |
+-------------+------------------------------------------------------------------------+
|<bc_id>      |The boundary flag identifier, an integer associated with the            |
|             |boundary location (node set in EXODUS II) in the problem                |
|             |domain on which the quantity is desired.                                |
+-------------+------------------------------------------------------------------------+
|<blk_id>     |An integer that designates the mesh block (material) from               |
|             |which the variable value should be taken. This has                      |
|             |implications for discontinuous variables on internal boundaries.        |
+-------------+------------------------------------------------------------------------+
|<species_id> |An integer that identifies the species number if a species              |
|             |variable is requested.                                                  |
+-------------+------------------------------------------------------------------------+
|<file_name>  |A character string corresponding to a file name into which              |
|             |the data should be printed.                                             |
+-------------+------------------------------------------------------------------------+

------------
**Examples**
------------

A simple example of using this card in context is shown below.
::

   Post Processing Data =
   DATA = VELOCITY2 1 1 0 data.out
   END OF DATA

-------------------------
**Technical Discussion**
-------------------------

If a fixed mesh or a subparametric mesh problem is being solved, the point coordinates
printed to the file will be the undeformed coordinates.



--------------
**References**
--------------

No References.