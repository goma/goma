3.2.15 Maximum path step
------------------------

**Maximum path step** = <float>

**Description/Usage**

This card is required for all continuation problems. <float> is the absolute value of the maximum allowable increment in the continuation parameter between steps.

**Examples**

To allow the continuation parameter step size to increase but not allow a step size larger than 2, use:

::

    Maximum path step = 2.0

**Technical Discussion**

For LOCA continuation problems, the function simple_step_control determines a step size for each step after the first. For continuation without LOCA, the function path_step_control performs this task. In either case, if a step size is computed which exceeds this limit, it is truncated to this value. For LOCA arc length continuation, an approximate limit is placed on the arc length step size.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
