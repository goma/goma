3.2.12 Maximum number of path steps
-----------------------------------

**Maximum number of path steps** = <integer>

**Description/Usage**

This card is required for all continuation problems. <integer> is the maximum number of continuation steps which may be taken. Continuation will stop after this number of steps even if the final parameter value is not reached.

**Examples**

To continue in a chosen continuation parameter from 1 to 5 at increments of 0.5 for a problem which has difficulty converging, and allow a maximum of 25 steps, use:

::

    Maximum number of path steps = 25

**Technical Discussion**

The initial continuation step, and every subsequent step attempt (successful or not), count toward this maximum.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
