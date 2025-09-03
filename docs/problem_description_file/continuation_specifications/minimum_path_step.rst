Minimum path step
------------------------

**Minimum path step** = <float>

**Description/Usage**

This card is required for all continuation problems. <float> is the absolute value of the smallest allowable increment in the continuation parameter between steps.

**Examples**

To set a minimum parameter step size of 0.1, use:

::

    Minimum path step = 0.01

This can also be entered in scientific notation as follows:

::

    Minimum path step = 1.0e-2

**Technical Discussion**

When convergence failure occurs on a continuation step, it is re-attempted at half of the previous step size until convergence is achieved. The value set for this cards places a limit on how small a step can be taken, and aborts continuation if the step size falls below this value. This can be used to avoid a large number of unnecessary step attempts, such as when the parameter reaches a problem stability limit.

When using arc length continuation, an approximate limit is placed on the arc length step size.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
