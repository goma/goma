3.2.5 Boundary condition data float tag
---------------------------------------

**Boundary condition data float tag** = <integer>

**Description/Usage**

This card is required for continuation problems when the specified Continuation Type is BC or AC. It is used to identify which of the floating-point inputs on the targeted BC (or user-supplied AC) card is to be used as the continuation parameter.

<integer>
    n, the zero-based numerical position (from left to right) of this value among the float inputs defined for this BC or user-supplied AC.
    
    -1, the constant value (e.g. volume, flux) is specified

**Examples**

A typical BC card may have a defined input format as follows:

::

    BC = {BC_name} <integer1> <integer2> <float1> <float2> <float3>

If the desired continuation parameter is <float3> above, then use:

::

    Boundary condition data float tag = 2

If you are using an augmenting condition and wish to continue in its constant value, use:

::

    Boundary condition data float tag = -1

If you have specified a user-defined AC in file user_ac.c which has a list of two floats and wish to continue in the first of these floats, use:

::

    Boundary condition data float tag = 0

**Technical Discussion**

There are many different input formats defined for BC cards. If in doubt, look up the correct format for the specific BC of interest.

Note that the BC ID cards for continuation are being "borrowed" for similar use by augmenting conditions, in lieu of creating a second set of cards.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
