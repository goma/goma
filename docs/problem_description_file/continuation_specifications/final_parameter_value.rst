Final parameter value
----------------------------

**Final parameter value** = <float>

**Description/Usage**

This card is required for all continuation problems. <float> is the parameter value at which continuation is to end.

**Examples**

To continue in a chosen parameter from 1 to 5 in increments of 0.5, use:

::

    Final parameter value = 5.0

**Technical Discussion**

The final parameter value may be higher or lower than the initial parameter value, allowing for continuation in either direction. The size of the final step will be adjusted if necessary to reach this exact value. When doing LOCA arc length continuation, however, this can only approximately be done.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
