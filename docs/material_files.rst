===============================
Material Files
===============================

The material (“mat”) file for *Goma* contains a description or specification of all the properties
required for the multi-physics capabilities of *Goma*. A separate *.mat* file must be developed for
each material present in each simulation. The *mat* file (see Figure 5) is split into seven sections:
(1) Physical Properties (Section 5.1), (2) Mechanical Properties and Constitutive Equations
(Section 5.2), (3) Thermal Properties (Section 5.3), (4) Electrical Properties (Section 5.4), (5)
Microstructure Properties (Section 5.5), (6) Species Properties (Section 5.6), and (7) Source
Terms (Section 5.7).

Each section in this chapter discusses a separate part of the material file specification and it
indicates the data cards or input records that may be used, followed by the options available for
each individual record (or line in the file) and the necessary input data/parameters. All input data
are specified in a free field format with successive data items separated by blanks or tabs. In this
version of the user’s manual, a new format has been instituted in which each record is presented in
a template structure. This template has eight parts: 1) a title, which is also the card name, 2) a
syntax, which is enclosed in a framed box and shows the proper contents of the card, 3) a
Description/Usage section, which presents the user options and descriptions of proper input
records, 4) an Example, 5) a Technical Discussion to provide relevant information to help the user
understand how to select from among various options or how to properly determine the desired
parameters, 6) a Theory to provide an understanding of the physics and mechanics that have been
implemented or are being exercised, 7) a FAQs section to present important user experience, and
8) a Reference section to identify citations and/or provide background information to the user.
This is a more lengthly but a more complete form for documenting and instructing users of *Goma*.

The syntax entry denotes a unique string for each input record which *Goma* parses in the input
file. All words in these unique strings are separated by a single white space and because the code
parses for these exact strings, the parser becomes case sensitive. The identifying string for a
particular specification is followed by an ‘=’ character. Following this character will be all
additional data for that record, if any. In the syntax box, this additional data is symbolically
represented by one or more variables with some appropriate delimiters. Typically, the user will
find a variable called *model_name* enclosed in **curly braces** ‘{}’; this would then be followed by
a description of specific options for *model_name* in the Description/Usage section. The curly
braces indicate a required input and that the user must select one of the offered options for
*model_name*. **Required parameters**, if any, for the model option are enclosed in **angle brackets** ‘<
>’, while **optional parameters** for model_name are enclosed in **square brackets** ‘[ ]’. Following
the ‘=’ character, the user may use white space freely between and among the remaining
parameters on the command line.

Figure 4 illustrates a typical material file. The section *headers*, e.g., “--- Physical Properties”, are
user comments that are not processed by the input parser. In all sections of this chapter,
*model_name* is a character string and *floating_point_const_list* is a list of floating point numbers
of arbitrary length separated by a comma or one or more white spaces. The remainder of this
chapter covers each card (line) of the material-description file in detail. For each parameter that is
not dimensionless, base units are indicated in square brackets ( [ ] ) at the end of the syntax line;
the base units are those indicated in the Nomenclature section of this document. Empty brackets (
[ ] ) denote dimensionless parameters, while those without units or brackets are simply model
names, other strings, or integers. Several model parameters, e.g., Diffusivity, where the model
options include other than the **CONSTANT** type with a single input value, identify the units as
[**varied**]. In these cases, the parameter units will be listed for the **CONSTANT** model option and
the units for individual input parameters will be identified in the parameter description.


.. figure:: /figures/340_goma_physics.png
	:align: center
	:width: 90%

	Sample material-description file format. Lines highlighted in bold-face type are required.

All property models will eventually have a **USER** and a **USER_GEN** option. When the former is
selected, the user must add the user model to the appropriate routine in the file *user_mp.c.* This
file contains a template to simplify the implementation of a model in a full-Newton context, but
has the restriction that none of the models can contain a dependence on gradients of variables. For
more complex models, which contain such dependencies, the user must resort to the more
sophisticated mechanism that comprise the routines in *user_mp_gen.c*

A relatively new capability/model available on many of the properties is a table-lookup feature.
That is, if the model is of type **TABLE**, then a linear or bilinear interpolation is used to extract the
material property value from a table of numbers representing the dependence. The best way to
explain this is with an example. Often times a property is dependent on temperature, or related
dependent variable. If discrete data is available of the property value at various temperatures, as
from a spreadsheet, then such a table can be read and with appropriate interpolation operations the
property value is determined. Throughout the material property options, the reader might see aat
**TABLE** option. The syntax for the input of that option is as follows:

<**Property name**> = **TABLE** *<integer1> <character_string1> [character_string2] {LINEAR |
BILINEAR} [integer2] [FILE = filenm]*

Here, the integers, character strings and floats are defined as follows:

<integer1> - the number N of columns in the table. The first N-1 columns are the values of the
independent variables (e.g. temperature, concentration, etc.) and the final Nth column is the
property value. This number is usually 2.

<character_string1> - Required variable name for first column. Valid variable names are
**TEMPERATURE, MASS_FRACTION, SPECIES, CAP_PRES, FAUX_PLASTIC,** and
**LOWER_DISTANCE**. The last three are specific to the Saturation model of porous flow,
the LAME Mu model, and the Lubrication Height function model, respectively.
Temperature and mass fraction dependence are available in all properties with a **TABLE** option
which make sense.

[character_string2] - Optional second variable name for bi-linear lookup dependence. This is
exploratory.

{*LINEAR | BILINEAR*} - type of interpolation

[*integer2*] - species number required only for MASS_FRACTION, SPECIES, and
FAUX_PLASTICITY variables.

[*FILE = <filenm>*] - The optional keyword ‘FILE=’ indicates that the table data is to be read
from a separate file identified by <filenm>. Each row of the table corresponds to one variable
value, and is input in free form CSV or space separated values. Note that if this ‘FILE=’ option is
not present then the data will be read from the input material file itself following the TABLE
model card. The end of the table is signaled by the keyword “END TABLE” (see example below).

Some examples are in order:

::

   Lame MU = TABLE 2 FAUX_PLASTIC 0 LINEAR FILE=stress_strain_comp.txt
   ...
   Lame MU = TABLE 2 TEMPERATURE LINEAR
   1. 293
   2. 300
   3. 425.
   END TABLE


Finally, before we get started, the following is an option added to allow existing Chemkin
material property databases to be read in, basically obviating the need to even read the material
(mat) file. The detailed description of input records provided in this chapter thus applies to the
case when the Default Database is set to GOMA_MAT.

.. toctree::
   :maxdepth: 1

   material_file/physical_properties
   material_file/mechanical_and_constitutive
   material_file/thermal_properties
   material_file/electrical
   material_file/microstructure
   material_file/species
   material_file/source_terms
   material_file/shell_equation
   material_file/moments


