3 Automatic Continuation
========================

Automatic continuation refers to the family of algorithms that allow tracking steady-state solution paths as a set of one or more parameters are varied. Goma is currently capable of zero order, first order, and arc length continuation in a single parameter set, and bifurcation tracking (turning point, pitchfork, or Hopf) in two parameter sets. Continuation can be carried out in four types of parameters:

**Material Properties**
    Material Property Tags are obtained from Goma source file mm_mp_const.h. A partial list is given in the previous chapter.

**Boundary Condition Floats**
    Continuation can be carried out in any boundary condition data float, i.e., the boundary condition parameters, making continuation in any boundary data and geometry simple. The definition of a boundary condition has the form:

    BC = BC_NAME BC_TYPE BC_ID DATAFLOAT1 DATAFLOAT2 . . .

    An example of continuation in a boundary condition: the effect of flow rate on a flow field can be tracked by continuing in the data float that represents velocity in a Dirichlet boundary condition.

**Material Property User Model Floats**
    This option is used when the user supplies a model for a material property in the file "user_mp.c" which includes a list of one or more input data floats, and one of these floats is the desired continuation parameter.

**Augmenting Conditions**
    There are two types of augmenting condition values which can be used as continuation parameters: the constant value (e.g. flux, volume) which is to be imposed, and one of a list of optional floats which may be specified with a user-defined augmenting condition (file "user_ac.c").

Following is a summary of solution prediction algorithms used in Goma's continuation schemes.

+------------------+------------------------------------------------------------------+
| Continuation     | Characteristic Algorithm                                         |
| Method           |                                                                  |
+==================+==================================================================+
| zero             | Single-parameter zeroth order continuation                       |
|                  |                                                                  |
|                  | x\ :sup:`PREDICTED`\ (λ\ :sub:`NEW`\ ) = x\ :sup:`OLD`\ (λ\     |
|                  | :sub:`OLD`\ )                                                    |
+------------------+------------------------------------------------------------------+
| first            | Single-parameter first order continuation                        |
|                  |                                                                  |
|                  | x\ :sup:`PREDICTED`\ (λ\ :sub:`NEW`\ ) = x\ :sup:`OLD`\ (λ\     |
|                  | :sub:`OLD`\ ) + Δλ\ :sub:`s`\ ∂x/∂λ                             |
+------------------+------------------------------------------------------------------+
| hzero            | Multi-parameter zeroth order continuation                        |
|                  |                                                                  |
|                  | x\ :sup:`PREDICTED`\ (λ\ :sub:`NEW`\ ) = x\ :sup:`OLD`\ (λ\     |
|                  | :sub:`OLD`\ )                                                    |
+------------------+------------------------------------------------------------------+
| hfirst           | Multi-parameter first order continuation                         |
|                  |                                                                  |
|                  | x\ :sup:`PREDICTED`\ (λ\ :sub:`NEW`\ ) = x\ :sup:`OLD`\ (λ\     |
|                  | :sub:`OLD`\ ) + Σ Δλ\ :sub:`j`\ ∂x/∂λ\ :sub:`j`\               |
+------------------+------------------------------------------------------------------+

These algorithms are available with or without the Library of Continuation Algorithms (LOCA), which also offers algorithms for arc length continuation and bifurcation tracking, as described above. Details of these algorithms can be found in the LOCA 1.0 manual (SAND 2002-0396).

3.1 Required Specifications in the Goma Input File
===================================================

A new section must be added to the Goma input file to identify the continuation method, continuation type, the continuation parameter, and other needed data. It is required only when automatic continuation is used in the numerical model and has the following form:

::

    ------------------------------------------------------------
    Continuation Specifications
    ------------------------------------------------------------
    Continuation = zero
    Continuation Type = MT
    Boundary condition ID = 4
    Boundary condition data float tag = 0
    Material id = 1
    Material property tag = 1700
    Material property tag subindex = 0
    Initial parameter value = 0.0
    Final parameter value = 800.0
    delta_s = 20.0
    Maximum number of path steps = 20
    Minimum path step = 1.0e-05
    Maximum path step = 100.0
    Continuation Printing Frequency = 1

The input records above control the continuation process. Specifically, this example indicates that zeroth-order continuation in material property density (tag 1700) will be simulated with Goma. (Boundary condition records must be present even though their input is ignored, and vice versa for material property records when the Continuation Type is BC.) The material subindex is not currently used. Density will be varied between 0 and 800 in material 1 using an initial path length/step size of 20.0 and minimum and maximum values of 1.0e-5 and 100.0, respectively. A maximum of 20 steps will be taken and every solution step will be written to the database. Each entry of the Continuation Specifications is discussed further in Section 3.2.

There are additional input cards which are either optional or required only in certain cases; these will also be described in Section 1.2.

.. include:: continuation_specifications/index.rst

3.3 Single Parameter Continuation via the Command Line
=======================================================

After a continuation run, Goma writes a text file called goma-cl.txt that lists the command line arguments for continuing directly from the command line. This makes it easy for writing scripts that can do multiple continuation sequences without having to provide multiple input files. There may be some benefit for use with other Goma wrap-arounds such as Dakota.

Continuation command line arguments take precedence over the values provided in the input file. Continuation is enabled even if Continuation = none is in the input file. For the Continuation Specifications section above (Section 3.2), goma-cl.txt contains the following after running the continuation sequence once:

::

    goma -a -i input -cb 0.000000e+00 -ce 8.000000e+02 \
    -cd 2.000000e+01 -cn 20 \
    -cm 0 -ct 2 \
    -c_mn 1 -c_mp 1700

These are the command line arguments required to replicate the main continuation commands specified in the input file. The backslashes indicate that all terms should appear on one line. All the continuation command line arguments (list by typing goma –help) are as follows:

::

    -cb FLT    Continuation: Start value
    -ce FLT    Continuation: Final value
    -cd FLT    Continuation: Path step, ds
    -cn INT    Continuation: Max number of path steps
    -cm INT    Continuation: Method
    -ct INT    Continuation: Type
    -c_bc INT  Continuation: Boundary condition ID
    -c_df INT  Continuation: BC Data Float ID
    -c_mn INT  Continuation: Material ID
    -c_mp INT  Continuation: Material property ID

The method and type flags are as follows: 0 for the method (cm) indicates zeroth order continuation and 1 first order continuation, while for the type, 1 indicates a boundary condition and 2 a material property. If the continuation parameter is a boundary condition data float, the Goma boundary condition and data float indices must be provided on the command line or in the input file. Similarly, if the continuation parameter is a material property, the material number and property tag must be provided on the command line or in the input file.

3.4 Multi-parameter Continuation with Hunting
==============================================

If the continuation method card (Section 3.2.1) is set to hzero or hfirst, that is

::

    Continuation = hzero/hfirst

the hunting continuation capability is used. This allows multi-parameter continuation in either boundary condition data floats or material properties or both. hzero and hfirst refer to the zeroth order and the first order methods, respectively. Both methods are enabled with a ramping feature so that all parameters can be ramped from a start to an end value in the specified number of steps. Hunting continuation can be used to continue in a single parameter instead of the standard zero and first methods. The maximum number of steps and printing frequency are specified in the Continuation Specifications section (3.2.12 and 3.2.18) described above.

When the hunting capability is being used, Goma ignores the Continuation Type, Boundary condition ID, Boundary condition data float tag, Material id, Material property tag, Material property tag subindex, Initial parameter value, Final parameter value, and delta_s cards in the Continuation Specifications section. (Those input records not ignored are the Continuation, Maximum number of path steps, Minimum and Maximum path steps, and Continuation Printing Frequency.)

A new section must be added to the Goma input file to identify the hunting continuation parameters, their start and end value, and step size information. It has the following form:

::

    ------------------------------------------------------------
    Hunting Specifications
    ------------------------------------------------------------
    Number of hunting conditions = -1
    HC = BC 4 0 1 0.0 1.0 1.0 0.0001 1000.0
    HC = MT 1 1700 1 0.0 100.0 10.0 0.0001 1000.0
    END OF HC

In the above example, two hunting conditions are specified. One indicates continuation in a boundary condition data float and another in the density. Other entries are described below.

3.5 Hunting Specifications
===========================

3.5.1 Number of hunting conditions
-----------------------------------

**Number of hunting conditions** = <int>

**Description/Usage**

This card is required when multiple parameter continuation is to be performed in either of these two cases: (1) The Goma hunting routine is invoked by setting the Continuation card to hzero or hfirst, or (2) Continuation is to be performed in LOCA (Continuation = loca) and the Number of continuation conditions card is set to -2, indicating that the continuation conditions are to be taken from HC cards instead of CC cards. In either case, if <int> is positive then <int> cards will be read below this card, and if <int> is -1 then the HC cards between this card and the END OF HC card will be counted and read.

**Examples**

Consider a continuation problem in which the continuation parameter is used by three different boundary conditions, hence three continuation or hunting conditions are required. For either of the two cases described above, use:

::

    Number of hunting conditions = 3

to be followed by three HC cards (one for each of these BC's).

Alternatively, count the cards by using:

::

    Number of hunting conditions = -1

to be followed by the same three HC cards, then an "END OF HC" card.

**Technical Discussion**

A backward compatibility routine has been provided such that LOCA can use either HC cards or CC cards. Accordingly, when HC cards are used (even with LOCA), the cards which identify and set values and step sizes for the continuation parameter in the Continuation Specifications section, while still required to be present, are overwritten with the values provided in the first HC card.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.

3.5.2 HC
---------

**HC** = {string} <int1> <int2> <int3> <float1> <float2> <float3> <float4> <float5>

**Description/Usage**

This card is required for each of two or more quantities which must be updated at each continuation step when the Goma hunting routine is invoked (Continuation = hzero or hfirst) or when LOCA is invoked (Continuation = loca) and Number of continuation conditions is set to -2 to indicate that HC cards (rather than CC cards) are to be read. The arguments contain the information required to uniquely identify the quantity to be updated and set its value range, initial step size, and step size limits. These arguments are as follows:

{string} indicates the type of quantity to be updated. The valid options are:

**BC**
    Boundary condition float argument.

**MT**
    Constant material property.

**AC**
    Augmenting condition (constant value or float argument).

**UM**
    User-defined material property model float argument.

NOTE: When "UM" is chosen, a fourth <int> argument will be needed!

The required number and meaning of the <int> entries depend on {string} as follows:

{string} = BC: Three integers
    int1    BCID - Zero-based position of relevant BC card.
    int2    DFID - Zero-based float argument number on this BC card.
    int3    Step control flag (see below).

{string} = MT: Three integers
    int1    MTID - One-based material index number.
    int2    MPID - Property tag number (assigned in "mm_mp_const.h").
    int3    Step control flag (see below).

{string} = AC: Three integers:
    int1    BCID - Zero-based position of relevant AC card.
    int2    DFID - Either zero-based data float (for user-supplied AC's) or -1 to indicate the AC constant value (e.g. volume, flux).
    int3    Step control flag (see below).

{string} = UM: Four integers
    int1    MTID - One-based material index number
    int2    MPID - Property tag number (assigned in "mm_mp_const.h").
    int3    MDID - Zero-based user model float argument number.
    int4    Step control flag (see below).

The step control flag is always the last integer entry on the card. Valid options are:

0   No step control
1   Use step control

When step control is used, the step size is set to a constant value equal to the range (end_value - start_value) divided by one less than the maximum number of steps. Otherwise, the step size is recalculated at each step based on the Newton convergence rate.

The float entries are as follows:

float1  start_value: Value on first continuation step.
float2  end_value: Value at the end of the continuation run.
float3  start_step: Step size taken after first continuation step.
float4  min_step_value: Minimum allowable step size.
float5  max_step_value: Maximum allowable step size.

If the step size calculated for the next step exceeds max_step_value, it is reset to this value. If it falls below min_step_value, continuation is aborted.

There are no defaults for any entries on this card.

**Examples**

Consider a continuation problem in which three sides of a box are to be held at a constant temperature T, which is to be increased from 50 to 100 in constant steps of 10. This is specified with Dirichlet BC cards at the top of the BC list.

The corresponding HC cards would then be used:

::

    HC = BC 0 0 1 50.0 100.0 10.0 10.0 10.0
    HC = BC 1 0 1 50.0 100.0 10.0 10.0 10.0
    HC = BC 2 0 1 50.0 100.0 10.0 10.0 10.0

To allow the step size to start at 10 and range from 5 to 20, use for the first card:

::

    HC = BC 0 0 0 50.0 100.0 10.0 5.0 20.0

and make the corresponding changes to the other two cards.

**Technical Discussion**

HC cards differ from CC cards in that step control is done independently for each card (thus care must be taken to ensure that the updated quantities change in consistent linear proportions to each other) and in that they can be used with or without LOCA.

An alternative to HC cards is to use user-defined continuation functions, which allow nonlinear functions to be specified (see file "user_continuation.c").

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.

3.5.3 END OF HC
---------------

**END OF HC**

**Description/Usage**

This card is the companion card to the Number of hunting conditions card. The END OF HC card signals the end of the Hunting Conditions Specifications section of the Goma input. When the Number of hunting conditions (= -1) is set to negative one, Goma will read all the hunting condition cards until this card is encountered in the input file. This card may omitted if the integer N on the Number of hunting conditions card is not -1.

**Examples**

Typical usage of this card is illustrated below:

::

    Number of hunting conditions = -1
    .
    .
    .
    END OF HC

**Technical Discussion**

See companion card Number of hunting conditions.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.

3.6 Continuation condition (CC) cards
=====================================

3.6.1 Purpose
--------------

Some continuation problems involve the requirement to update a number of different conditions of the problem each time a continuation step is taken. For example, if there are three boundary condition cards which have a common float input, then each must be updated at every step. It may also be necessary to match some quantity on both sides of a material interface. Or, it may be of interest to examine the response of a system to two variables changing simultaneously, even if the physics of the problem does not demand it.

Hunting conditions, as discussed in Section 3.5, provide a mechanism for fulfilling this need by providing ID and step information for each quantity to be updated (using HC cards). However, this requires the use of one of two designated Goma hunting algorithms (hzero and hfirst), which are handled separately from all other continuation algorithms. In order to have the ability to handle multiple updates per step with the more advanced algorithms in LOCA, it is necessary to integrate this feature with each of the algorithms. For many such problems, the quantities to be updated maintain a linear relationship to each other (i.e. y=mx+b) as they change, either by design or by necessity. In fact, this linear relationship is required for the hunting algorithms to work.

3.6.2 Usage
------------

The approach taken for LOCA is to designate one of the update quantities as the "master" and place the ID and step information for it in the cards pertaining to a single continuation parameter (which the hunting algorithms require but do not use). Any other update quantities, or "continuation conditions", are considered "slaves" and identified by inserting a CC card, which need only include the necessary ID information and a means to specify the linear relationship between changes in it and changes in the "master", as described in Section 3.2. Thus, there is one fewer CC card required. A global variable indicates how many continuation conditions are in effect, so that when an update call comes from LOCA, the correct number of parameter update calls are made (even if only one).

Here is how the previous example continuation section could be restated using CC cards:

::

    ------------------------------------------------------------
    Continuation Specifications
    ------------------------------------------------------------
    Boundary condition ID = 4
    Boundary condition data float tag= 0
    ...
    Initial parameter value= 0.0
    Final parameter value= 10.0
    delta_s= 1.0
    Minimum path step= 0.0001
    Maximum path step= 1000.0
    LOCA method= zero
    Number of continuation conditions= 2 (or -1)
    CC = MT 1 1700 2 0.0 10.0
    END OF CC

The CC card here is set up such that MT 1 ID 1700 starts at 0.0 and changes by 10x whenever BC 4 float 0 changes by x. Note also that the number of continuation conditions is 2, but only one CC card is required because the information for the first (the "master") is taken from the preceding cards.

To maintain backward compatibility, LOCA can use an existing set of hunting (HC) cards in lieu of converting them to CC cards, and will extract the information it needs properly. To do so, it is only necessary to include the card:

::

    Number of continuation conditions= -2

which will tell the input parser to read hunting cards even if a hunting algorithm is not specified by the Continuation card.

3.7 User-defined continuation conditions
========================================

3.7.1 Purpose
--------------

When each quantity which must be updated during continuation bears a linear (y=b+mx) relationship to the designated continuation parameter, then these relationships can be adequately specified with either HC or CC cards. The CC and TC cards can also be used to specify a two-term polynomial relationship of the type y=b+mx^n, or relationships involving the trigonometric functions (e.g. sin, cos) of an angular continuation parameter. However, some problems may require more complicated relationships between these values, such as multiple-term polynomials, exponentials, etc. In other cases, the previous solution may be required to obtain new values. Such requirements would then require a mechanism for the user to construct custom continuation conditions for the problem at hand. This mechanism is provided in the file user_continuation.c.

3.7.2 Usage
------------

The file user_continuation.c contains two template functions for the user to supply custom continuation conditions: update_user_parameter (for the primary continuation parameter) and update_user_TP_parameter (for the second parameter used in LOCA bifurcation tracking algorithms). Detailed instructions and example entries are provided in the comments of this file.

To use this feature, first identify all quantities q\ :sub:`i`\ which are to be updated at each continuation step. Then designate a single scalar continuation parameter λ such that each update quantity can be readily expressed as a function of λ (and possibly the solution), viz.

q\ :sub:`i`\ = f(λ, x, ẋ, x\ :sub:`AC`\ )                             (3-3)

where x is the current solution vector, xdot is its time derivative, and x_AC is the vector of extra unknowns when augmenting conditions are used. Note that λ may or may not be one of the update quantities q\ :sub:`i`\ . It may be used, for example, as a "progress parameter", as theta is for hunting, or it can be a dimensionless group (such as Reynolds number) which is based on the input quantities but may not be equal to any of them. When user-defined continuation is specified for a parameter, the parameter ID card (e.g. BCID, MTID) card inputs are not used (but must still be present). However, the step inputs (Initial and final values, starting, minimum, and maximum step size) are used for λ. All other necessary information will be placed in the relevant user continuation function and compiled into the code; no HC/CC/TC cards are needed.

The procedure for supplying the continuation conditions consists of three steps for each quantity: 1) providing the necessary type and ID information, 2) defining the new value in terms of λ and the solution, and 3) pasting in a standard update call. This may be done for either or both parameters for LOCA bifurcation tracking algorithms. Once these functions are entered, it is then necessary to recompile Goma to include them. User-defined continuation is then invoked by specifying continuation type "UF" (user functions) in the input file and providing the number of continuation functions (conditions) entered for the relevant parameter. For example, to continue in web speed for a slot coater, the relevant update quantities are the two Dirichlet BC's for the web and the fluid in contact with it, and since it is necessary to vary the inflow proportionally to maintain constant downstream film thickness, each of the three floats of a GD_PARAB BC also. Thus, there are five quantities. Here is how the functions would be entered when λ is the web speed:

::

    {
     static int first_cp = 1;
     int first_tp = -1;
     int n=0;
     int Type, BCID, DFID, MTID, MPID, MDID;
     double value;

    /* Declare any additional variables here */
     double f2, f3, f4;

    /* If using this function, comment out this line. */
    /*EH(-1, "No user continuation conditions entered!");*/

    /* Evaluate any intermediate quantities and/or assign constants here */
     f2 = -9.684754406;
     f3 = 318.0179944;
     f4 = -2486.458128;

    /* Enter each continuation condition in sequence in this space */

    /* Condition 0: Slip velocity */
    /* ID's */
     Type = BC;
     BCID = 13;
     DFID = 1;
    /* Value */
     value = lambda;
    /* Update call - copy from example */
     n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
     MTID, MPID, MDID, value, cx, exo, dpi);

    /* Condition 1: Fluid velocity along web */
    /* ID's */
     Type = BC;
     BCID = 14;
     DFID = 0;
    /* Value */
     value = lambda;
    /* Update call */
     n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
     MTID, MPID, MDID, value, cx, exo, dpi);

    /* Condition 2: Inlet velocity GD_PARAB first float */
    /* ID's */
     Type = BC;
     BCID = 11;
     DFID = 0;
    /* Value */
     value = f2 * lambda;
    /* Update call */
     n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
     MTID, MPID, MDID, value, cx, exo, dpi);

    /* Condition 3: Inlet velocity GD_PARAB second float */
    /* ID's */
     Type = BC;
     BCID = 11;
     DFID = 1;
    /* Value */
     value = f3 * lambda;
    /* Update call */
     n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
     MTID, MPID, MDID, value, cx, exo, dpi);

    /* Condition 4: Inlet velocity GD_PARAB third float */
    /* ID's */
     Type = BC;
     BCID = 11;
     DFID = 2;
    /* Value */
     value = f4 * lambda;
    /* Update call */
     n = do_user_update(n, first_cp, first_tp, Type, BCID, DFID,
     MTID, MPID, MDID, value, cx, exo, dpi);

    /* Done */
     first_cp = 0;
     return;
    }

Then, to invoke this function, include the following Continuation input cards:

::

    Continuation Type = UF
    Number of user continuation functions = 5