delta_s
--------------

**delta_s** = <float>

**Description/Usage**

This card is required for all continuation problems. It specifies the size of the first parameter step to be taken; viz. λ\ :sub:`1`\ = λ\ :sub:`0`\ + delta_s.

<float>
    the size of the first parameter step.

**Examples**

To continue in a chosen parameter from 1 to 5 at increments of 0.5, use:

::

    delta_s = 0.5

**Technical Discussion**

If a continuation step fails when using LOCA with a constant specified step, the step size is initially halved and then allowed to increase back up to delta_s over a number of steps.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
