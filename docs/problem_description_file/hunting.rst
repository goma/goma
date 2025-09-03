Hunting Specifications
##########################

The cards in this section are used to direct multiparameter continuation with the hunting
technique, which is a linear, multiparameter capability. The user is referred to discussions for the
*Continuation Specifications* for the important details of Continuation. As is true for the
*Continuation Specifications*, this entire section is completely optional. *Hunting Specification*
cards are used in conjunction with *Continuation Specifications.*

Input specifications for this section of input records is discussed in a separate, comprehensive
manual (Gates, et. al., 2000); an update to this manual has been completed during the summer of
2006 (Labreche, et. al., 2006).

Number of hunting conditions
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

