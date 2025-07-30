3.2.3 Number of user continuation functions
-------------------------------------------

**Number of user continuation functions** = <integer>

**Description/Usage**

This card is required for all continuation problems when the selected Continuation Type is "UF".

<integer>
    the number of functions provided by the user in the function "update_user_parameter" (in file user_continuation.c). These functions identify quantities (e.g. BC floats, material properties) which must be updated on each continuation step and then calculates new values for them on each parameter update.

**Examples**

If there are three quantities to be updated at each continuation parameter step, use:

::

    Number of user continuation functions = 3

**Technical Discussion**

In order to use the "UF" continuation type, a list of functions which calculate each update value as a function of the continuation parameter λ must be provided and compiled in file "user_continuation.c". This approach differs from using multiple continuation conditions (specified by CC cards) in that the continuation parameter specified in the input deck is not used directly, but instead through these functions, parameter updates are performed. Thus, the BC/MT ID cards are not actually used (as this information is provided along with the functions), but the range and step size input cards are used. λ can then be used as a progress parameter (i.e. range from 0 to 1) as long as the user-provided functions are based on the λ range specified in these cards.

The number of user continuation functions must be provided in advance so that the correct number of copies of the User_Continuation_Info structure (cpuc) can be allocated.

The file user_continuation.c provides templates for two different user continuation function sets: update_user_parameter (for the first parameter) and update_user_TP_parameter (for the second parameter). The latter is used only with LOCA bifurcation tracking algorithms; a different card is used to enumerate these functions.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
