**************
**Initialize**
**************

::

   Initialize = {char_string} <integer> <float> [varies]

-----------------------
**Description / Usage**
-----------------------

This optional card provides a mechanism to set one of the field variables to a constant
value within the current material block. Definitions of the input parameters are as
follows:

+--------------------------+-------------------------------------------------------------------------------------+
|<char_string>             |Permissible values for this input string are any variable names identified in source |
|                          |file rf_fem_const.h beginning at the section labeled Variable Names of unknowns,     |
|                          |though they should be active in this material block. Examples include, but are not   |
|                          |limited to, the following:                                                           |
|                          |                                                                                     |
|                          |**VELOCITY1, VELOCITY2, VELOCITY3 (V123),MESH_DISPLACEMENT (MD123),                  |
|                          |SOLID_DISPLACEMENT (SD123), MASS_FRACTION, TEMPERATURE, PRESSURE,VOLTAGE, FILL, LS,  |
|                          |POLYMER_STRESS (6 components, 8 modes), VELOCITY_GRADIENT (9 components), SHEAR_RATE,|
|                          |VOLF_PHASE (6 phases), POR_LIQ_PRES, POR_GAS_PRES, POR_POROSITY, POR_SATURATION,     |
|                          |POR_LAST, LAGR_MULT (LM123), SURF_CHARGE, EXT_VELOCITY, EFIELD(123), SHELL           |
|                          |(4 variables), SPECIES (7 variables).**                                              |
|                          |                                                                                     |
|                          |*Note: for a comprehensive list of initializable variables, consult Volume 1         |
|                          |“Initialize” card.*                                                                  |
+--------------------------+-------------------------------------------------------------------------------------+
|<integer>                 |Species number to be initialized if the value of {char_string} is one of the         |
|                          |**SPECIES** variables (see Technical Discussion); otherwise, set <integer> to zero.  |
+--------------------------+-------------------------------------------------------------------------------------+
|<float>                   |Value to which the variable should be initialized.                                   |
+--------------------------+-------------------------------------------------------------------------------------+

*Multiple applications of this card are valid; Goma automatically counts the number of Initialize cards.*

------------
**Examples**
------------

Following is a sample card:

::

   INITIALIZE = POLYMER_STRESS11 0 1.25E4

-------------------------
**Technical Discussion**
-------------------------

This card provides the means to set initial values for any of the field variables in the
element block for a particular material. Since the setting of variables initialized on this
card takes place after reading the initial guess (see function init_vec in file rf_util.c), it
can be used to override the value in the initial guess file.

In order to set a field to a specific value over the entire problem domain, a similar
*Initialize* capability is provided as a global variable in the *General Specifications*
section of the *Goma* input file. Please check in the Problem Description section of this
manual.

Note, the **SPECIES_UNK** variables are **NOT** used to initialize any of the species
variables. Rather, the special definition called **MASS_FRACTION**
is the variables used in Goma input or mat
files for this input record. Multiple species are initialized by combining 
**MASS_FRACTION** with the second parameter (<integer>) on this card. These cards are
particularly handy for mass transfer problems, where the initial conditions need to
specify different concentrations of the same species in different materials.

Note: for a comprehensive list of initializable variables, consult Volume 1 “Initialize”
card.



