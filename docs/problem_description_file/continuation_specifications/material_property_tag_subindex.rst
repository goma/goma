Material property tag subindex
------------------------------------

**Material property tag subindex** = <integer>

**Description/Usage**

This card is required for continuation problems when the specified Continuation Type is "UM". It identifies the continuation parameter by its zero-based position (left to right) within a list of float parameters specified for a user-defined property of a material.

**Examples**

Consider a case where the desired continuation parameter is the temperature dependence of surface tension. Currently, doing this in Goma requires a surface tension model to be supplied in "user_mp.c" such as:

::

    sigma = param[0] - param[1] * T;
    dsigma_dT = -param[1];

Here, the continuation parameter would be param[1]. To employ this model, the \*.mat file would include the card:

::

    Surface Tension = USER <float1> <float2>

To specify param[1] as the continuation parameter, use:

::

    Material property tag subindex = 1

**Technical Discussion**

The number of floats assigned for a given user property model is determined by counting the number supplied with the property model card in the \*.mat file. If N floats are given there, then they are assigned indices from 0 to N-1 in the param[] array. Accordingly, the subindex number must not exceed N-1, or an error will occur.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
