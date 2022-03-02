*********
Viscosity
*********

::

   Viscosity = {model_name} {float_list} [M/Lt]

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the viscosity model for the liquid constitutive equation (see
Sackinger et al., 1995). Definitions of the input parameters are as follows:

+-----------------+----------------------------------------------------------------------------------------------------+
|{model_name}     |The name of the viscosity model, which can be one of the following: **CONSTANT, USER, USER_GEN**, or|
|                 |**FILL, LEVEL_SET, CONST_PHASE_FUNCTION**.                                                          |
+-----------------+----------------------------------------------------------------------------------------------------+
|{float_list}.    |One or more floating point numbers (<float1> through <floatn>) whose values are determined by the   |
|                 |selection for {model_name}. These are identified in the discussion of each model below. Note that   |
|                 |not all models employ a {float_list}.                                                               |
+-----------------+----------------------------------------------------------------------------------------------------+

Thus,

+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**CONSTANT** <float1>                                                                  |This option specifies a constant viscosity for a Newtonian fluid. The      |
|                                                                                       |{float_list} has a single value:                                           |
|                                                                                       |                                                                           |
|                                                                                       | * <float1> - value of viscosity                                           |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**USER** <float1>... <floatn                                                           |This option specifies that the viscosity will be given by a user-defined   |
|                                                                                       |model; the model must be incorporated into *Goma* by modifying function    |
|                                                                                       |“usr_viscosity” in file user_mp.c. The model parameters are entered in the |
|                                                                                       |{float_list} as <float1> through <floatn> and passed to the routine as an  |
|                                                                                       |array.                                                                     |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**USER_GEN** <float1>... <floatn>                                                      |This option specifies that the viscosity will be given by a generalized    |
|                                                                                       |user-defined model. This user-defined model must be incorporated by        |
|                                                                                       |modifying the routine “usr_viscosity_gen” in the file user_mp_gen.c. Any   |
|                                                                                       |number of parameters can be passed (via <float1> through <floatn>) in here.|
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**FILL** <float1> <float2>                                                             |The {float_list} for this option requires two values. It invokes a FILL    |
|                                                                                       |dependent viscosity that is set to the value of *float1* if the FILL       |
|                                                                                       |variable is 1 and *float2* if the FILL variable is 0.                      |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**LEVEL_SET** <float1> <float2> <float3>                                               |This model is used to vary the viscosity in the flow region when a level   |
|                                                                                       |set function is used to track the boundary between two fluids using level  |
|                                                                                       |set interface tracking. This choice assures a smooth transition in density |
|                                                                                       |across the zero level set contour. The {float_list} contains three values  |
|                                                                                       |for this model, where:                                                     |
|                                                                                       |                                                                           |
|                                                                                       | * <float1> Fluid viscosity in the negative regions of the level set       |
|                                                                                       |   function, μ-                                                            |
|                                                                                       | * <float2> Fluid viscosity in the positive regions of the level set       |
|                                                                                       |   function, μ+                                                            |
|                                                                                       | * <float3> Length scale over which the transition occurs, :math:`\alpha` .|
|                                                                                       |   If this parameter is set to zero, it will default to one-half the Level |
|                                                                                       |   Set Length Scale value specified elsewhere in the input deck.           |
|                                                                                       |                                                                           |
|                                                                                       |**Note**: a better way to specify the identical viscosity model is to make |
|                                                                                       |use of the 2nd Level Set Viscosity card documented also in this manual.    |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**CONST_PHASE_FUNCTION** <floatlist> <float1> <float2>                                 |This model is used to vary the viscosity in the flow regime when phase     |
|                                                                                       |functions are used to track the motion of muliple phases. This choice      |
|                                                                                       |assures a smooth transition in viscosity across the phase boundaries. The  |
|                                                                                       |{float_list} contains a variable number of values that depend on the number| 
|                                                                                       |phase functions being tracked, where:                                      |
|                                                                                       |                                                                           |
|                                                                                       | * <floatlist> list of float variables equal to the number of phase        |
|                                                                                       |   functions. These are the constant viscosities associated with each      |
|                                                                                       |   phase in order from 1 to number of phase functions.                     |
|                                                                                       | * <float1> Length scale over which the transition between one phases      |
|                                                                                       |   viscosity value to theother occurs, :math:`\alpha` . If this parameter  |
|                                                                                       |   is set to zero, it will default to one-half the Level Set Length Scale  |
|                                                                                       |   value already specified.                                                |
|                                                                                       | * <float3> The “null” value for viscosity. This is the value for viscosity|
|                                                                                       |   which will be assigned to those regions of the flow where all the phase |
|                                                                                       |   functions are less than or equal to zero.                               |
|                                                                                       |                                                                           |
|                                                                                       |The user should examine the **CONST_PHASE_FUNCTION** density model for a   |
|                                                                                       |detailed description of the relations used to compute viscosity with this  |
|                                                                                       |model. That model refers to densities but the same equations apply if      |
|                                                                                       |viscosities are exchanged for densities.                                   |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+
|**TABLE** <integer1> <character_string1> {LINEAR | BILINEAR} [integer2] [FILE = filenm]|Please see discussion at the beginning of the material properties chapter 5|
|                                                                                       |for input description and options. Currently the only valid options for    |
|                                                                                       |character_string1 is **TEMPERATURE** and **MASS_FRACTION**.                |
+---------------------------------------------------------------------------------------+---------------------------------------------------------------------------+

------------
**Examples**
------------

The following is a sample card that sets the viscosity to USER:

::

   Viscosity = USER 1. 1. 1. 1. 1.

::

   Viscosity = LEVEL_SET 0.083 0.0001 0.1

-------------------------
**Technical Discussion**
-------------------------

The viscosity specified by this input card is used with the **NEWTONIAN** *Liquid
Constitutive Equation*.



--------------
**References**
--------------

Sackinger, P. A., Schunk, P. R. and Rao, R. R. 1995. "A Newton-Raphson Pseudo-Solid
Domain Mapping Technique for Free and Moving Boundary Problems: A Finite
Element Implementation", J. Comp. Phys., 125 (1996) 83-103.
