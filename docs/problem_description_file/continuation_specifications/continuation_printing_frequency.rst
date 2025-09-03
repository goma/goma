Continuation Printing Frequency
--------------------------------------

**Continuation Printing Frequency** = <integer>

**Description/Usage**

This card is required for all continuation problems. It is used to specify that the solution is to be written to the EXODUS (*.exoII) and ASCII (*.dat) output files after a specified number of steps.

<integer>
    N, frequency of writing a continuation step

**Examples**

To write the solution after every step, use:

::

    Continuation Printing Frequency = 1

To write the solution after every third step, use

::

    Continuation Printing Frequency = 3

To output the solution only after the last step, use:

::

    Continuation Printing Frequency = N+1

where N is the maximum number of path steps (the last step is always written).

**Technical Discussion**

This card can be used to avoid generating very large output files, or when only the final step is of interest.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
