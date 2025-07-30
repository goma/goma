3.2.2 Continuation Type
-----------------------

**Continuation Type** = {BC | MT | AC | UM | UF | AN}

**Description/Usage**

This card is required for all continuation methods (including LOCA methods) to identify the type of continuation parameter. Valid options are:

**BC**
    Boundary condition input float.

**MT**
    Pre-defined constant material property.

**AC**
    Augmenting condition (constant value or input float).

**UM**
    User-defined material property model input float.

**UF**
    User-defined continuation function list.

**AN**
    Angular continuation parameter

The option selected determines which of the subsequent cards in this section are necessary to uniquely identify the continuation parameter.

**Examples**

For a parameter which appears in a boundary condition card, use:

::

    Continuation Type = BC

For a constant material property, use:

::

    Continuation Type = MT

For a parameter which appears in an augmenting condition card, use:

::

    Continuation Type = AC

For a parameter used in a user-defined material property model, use:

::

    Continuation Type = UM

For a user-supplied continuation function set, use:

::

    Continuation Type = UF

If the continuation parameter will be an angle and trigonometric functions of it will be needed, use:

::

    Continuation Type = AN

**Technical Discussion**

To use the BC option, a BC card must be provided for the relevant boundary conditions on the Boundary Condition Specification section.

To use the MT option, a tag for the relevant property must be defined in the file mm_mp_const.h.

To use the AC option, an AC card must be provided for the relevant augmenting condition in the Augmenting Conditions Specifications section.

To use the UM option, a property tag must be defined (as for MT), a user model and parameter list must be specified for the relevant property in file user_mp.c and this USER model must be specified with the correct number of parameter values in the relevant material (*.mat) file.

To use the UF option, a list of user-defined continuation conditions (functions) must be provided in the function update_user_parameter, which is in the file user_continuation.c.

The AN (angular parameter) option works differently in that the quantities to be updated use trigonometric functions (e.g. sin, cos) of this angle, rather than the angle itself. These functions are specified in the CC or TC cards which follow - these cards have completely different interpretations than with any of the other Continuation Type options (see the entry for the CC card). At least one such card must be provided.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
