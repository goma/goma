**************
**Initialize**
**************

::

	Initialize = {char_string} <integer> <float>       [units vary]

-----------------------
**Description / Usage**
-----------------------

This optional card provides a mechanism to set one of the field variables to a constant
value across the *whole* domain. Definitions of the input parameters are as follows:

{char_string}
    Permissible values for this input string are any variable names identified
    in source file rf_fem_const.h beginning at the section labeled Variable
    Names of unknowns. Examples include, but are not limited to, the following
    (note the shorthand notation for components):

    ``VELOCITY1, VELOCITY2, VELOCITY3 (V123),``
    ``MESH_DISPLACEMENT (MD123), SOLID_DISPLACEMENT``
    ``(SD123), MASS_FRACTION, TEMPERATURE, PRESSURE,``
    ``VOLTAGE, FILL, LS, POLYMER_STRESS (6 components,``
    ``8 modes), VELOCITY_GRADIENT (9 components),``
    ``SHEAR_RATE, VOLF_PHASE (6 phases), POR_LIQ_PRES,``
    ``POR_GAS_PRES, POR_POROSITY, POR_SATURATION,``
    ``POR_LAST, LAGR_MULT (LM123), SURF_CHARGE,``
    ``EXT_VELOCITY, EFIELD(123), SHELL (4 variables),``
    ``SPECIES (7 variables).``

    *For a more comprehensive list, see Technical discussion below.*

<integer>
    Species number to be initialized if the value of {char_string} is one of
    the SPECIES variables (see Technical Discussion); otherwise, set <integer>
    to zero.

<float>
    Value to which the variable should be initialized.

Multiple applications of this card are valid; *Goma* automatically counts the number of
*Initialize* cards.

------------
**Examples**
------------

Following is a sample card:
::

	Initialize = VELOCITY1 0 0.

-------------------------
**Technical Discussion**
-------------------------

This card provides the means to globally set (i.e., the entire problem domain) initial
values for any of the field variables. Since the setting of variables initialized on this
card takes place after reading the initial guess (see function *init_vec* in file *rf_util.c*), it can be used to override the value in the *Initial Guess* file.

In order to set a field to a specific value in a particular material only, a similar *Initialize*
capability is provided within each material block. Please check in the Material Files
section of this manual.

Note, the SPECIES_UNK variables are **NOT** used to initialize any of the species
variables. Rather, the special definition called **MASS_FRACTION**
representing the various Species Types, is the variables used in Goma input or mat
files for this input record. Multiple species are initialized by combining one of these
variable types with the second parameter (<integer>) on this card.

The comprehensive list of keyword variable names can be found in *mm_input_util.c*, if
you have access to GOMA source code. Search for the function *variable_string_to_int*.
A snapshot of the initialize-able variables in that routine is shown here:

::

    var = VELOCITY1;
    var = VELOCITY2;
    var = VELOCITY3;
    var = TEMPERATURE;
    var = MASS_FRACTION;
    var = MESH_DISPLACEMENT1;
    var = MESH_DISPLACEMENT2;
    var = MESH_DISPLACEMENT3;
    var = PRESSURE;
    var = POLYMER_STRESS11;
    var = POLYMER_STRESS12;
    var = POLYMER_STRESS13;
    var = POLYMER_STRESS22;
    var = POLYMER_STRESS23;
    var = POLYMER_STRESS33;
    var = SOLID_DISPLACEMENT1;
    var = SOLID_DISPLACEMENT2;
    var = SOLID_DISPLACEMENT3;
    var = VELOCITY_GRADIENT11;
    var = VELOCITY_GRADIENT12;
    var = VELOCITY_GRADIENT13;
    var = VELOCITY_GRADIENT21;
    var = VELOCITY_GRADIENT22;
    var = VELOCITY_GRADIENT23;
    var = VELOCITY_GRADIENT31;
    var = VELOCITY_GRADIENT32;
    var = VELOCITY_GRADIENT33;
    var = VOLTAGE;
    var = FILL;
    var = SHEAR_RATE;
    var = PVELOCITY1;
    var = PVELOCITY2;
    var = PVELOCITY3;
    var = POLYMER_STRESS11_1;
    var = POLYMER_STRESS12_1;
    var = POLYMER_STRESS22_1;
    var = POLYMER_STRESS13_1;
    var = POLYMER_STRESS23_1;
    var = POLYMER_STRESS33_1;
    var = POLYMER_STRESS11_2;
    var = POLYMER_STRESS12_2;
    var = POLYMER_STRESS22_2;
    var = POLYMER_STRESS13_2;
    var = POLYMER_STRESS23_2;
    var = POLYMER_STRESS33_2;
    var = POLYMER_STRESS11_3;
    var = POLYMER_STRESS12_3;
    var = POLYMER_STRESS22_3;
    var = POLYMER_STRESS13_3;
    var = POLYMER_STRESS23_3;
    var = POLYMER_STRESS33_3;
    var = POLYMER_STRESS11_4;
    var = POLYMER_STRESS12_4;
    var = POLYMER_STRESS22_4;
    var = POLYMER_STRESS13_4;
    var = POLYMER_STRESS23_4;
    var = POLYMER_STRESS33_4;
    var = POLYMER_STRESS11_5;
    var = POLYMER_STRESS12_5;
    var = POLYMER_STRESS22_5;
    var = POLYMER_STRESS13_5;
    var = POLYMER_STRESS23_5;
    var = POLYMER_STRESS33_5;
    var = POLYMER_STRESS11_6;
    var = POLYMER_STRESS12_6;
    var = POLYMER_STRESS22_6;
    var = POLYMER_STRESS13_6;
    var = POLYMER_STRESS23_6;
    var = POLYMER_STRESS33_6;
    var = POLYMER_STRESS11_7;
    var = POLYMER_STRESS12_7;
    var = POLYMER_STRESS22_7;
    var = POLYMER_STRESS13_7;
    var = POLYMER_STRESS23_7;
    var = POLYMER_STRESS33_7;
    var = SPECIES_MASS_FRACTION;
    var = SPECIES_MOLE_FRACTION;
    var = SPECIES_VOL_FRACTION;
    var = SPECIES_DENSITY;
    var = SPECIES_CONCENTRATION;
    var = SPECIES_CAP_PRESSURE;
    var = SPECIES_UNDEFINED_FORM;
    var = POR_LIQ_PRES;
    var = POR_GAS_PRES;
    var = POR_POROSITY;
    var = POR_TEMP;
    var = POR_SATURATION;
    var = VORT_DIR1;
    var = VORT_DIR2;
    var = VORT_DIR3;
    var = CURVATURE;
    var = BOND_EVOLUTION;
    var = SURF_CHARGE;
    var = EXT_VELOCITY;
    var = EFIELD1;
    var = EFIELD2;
    var = EFIELD3;
    var = ENORM;
    var = NORMAL1;
    var = NORMAL2;
    var = NORMAL3;
    var = SHELL_CURVATURE;
    var = SHELL_TENSION;
    var = SHELL_X;
    var = SHELL_Y;
    var = SHELL_USER;
    var = PHASE1;
    var = PHASE2;
    var = PHASE3;
    var = PHASE4;
    var = PHASE5;
    var = SHELL_ANGLE1;
    var = SHELL_ANGLE2;
    var = SHELL_SURF_DIV_V;
    var = SHELL_SURF_CURV;
    var = N_DOT_CURL_V;
    var = GRAD_V_DOT_N1;
    var = GRAD_V_DOT_N2;
    var = GRAD_V_DOT_N3;
    var = ACOUS_PREAL;
    var = ACOUS_PIMAG;
    var = ACOUS_ENERGY;
    var = POR_SINK_MASS;
    var = VORT_DIR1
    var = VORT_DIR2
    var = VORT_DIR3
    var = VORT_LAMBDA
    var = CURVATURE
    var = LAGR_MULT1
    var = LAGR_MULT2
    var = LAGR_MULT3
    var = BOND_EVOLUTION
    var = SURF_CHARGE
    var = EXT_VELOCITY
    var = EFIELD1
    var = EFIELD2
    var = EFIELD3
    var = ENORM
    var = NORMAL1
    var = NORMAL2
    var = NORMAL3
    var = SHELL_CURVATURE
    var = SHELL_TENSION
    var = SHELL_X
    var = SHELL_Y
    var = SHELL_USER
    var = PHASE1
    var = PHASE2
    var = PHASE3
    var = PHASE4
    var = PHASE5
    var = SHELL_ANGLE1
    var = SHELL_ANGLE2
    var = SHELL_SURF_DIV_V
    var = SHELL_SURF_CURV
    var = N_DOT_CURL_V
    var = GRAD_S_V_DOT_N1
    var = GRAD_S_V_DOT_N2
    var = GRAD_S_V_DOT_N3
    var = ACOUS_PREAL
    var = ACOUS_PIMAG
    var = SHELL_DIFF_FLUX
    var = SHELL_DIFF_CURVATURE
    var = SHELL_NORMAL1
    var = SHELL_NORMAL2
    var = ACOUS_REYN_STRESS
    var = SHELL_BDYVELO
    var = SHELL_LUBP
    var = LUBP
    var = SHELL_FILMP
    var = SHELL_FILMH
    var = SHELL_PARTC
    var = SHELL_SAT_CLOSED
    var = SHELL_PRESS_OPEN
    var = SHELL_TEMPERATURE
    var = SHELL_DELTAH
    var = SHELL_LUB_CURV
    var = SHELL_SAT_GASN
    var = SHELL_SHEAR_TOP
    var = SHELL_SHEAR_BOT
    var = SHELL_CROSS_SHEAR
    var = MAX_STRAIN
    var = CUR_STRAIN
    var = LUBP_2
    var = SHELL_PRESS_OPEN_2
    var = SHELL_LUB_CURV_2



