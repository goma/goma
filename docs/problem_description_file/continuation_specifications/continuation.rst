Continuation
-------------------

**Continuation** = {zero | first | hzero | hfirst | loca}

**Description/Usage**

This card is required only for continuation problems. It is used to specify the continuation method to be used. Valid options are:

**zero**
    Zero-order continuation.

**first**
    First-order continuation.

**hzero**
    Zero-order hunting.

**hfirst**
    First-order hunting.

**loca**
    Use the Library of Continuation Algorithms (LOCA).

If this card is not present, then none of the other Continuation Specifications cards in this section or Hunting Specifications cards in the Hunting section are needed.

**Examples**

Sample card for zero-order single-parameter continuation:

::

    Continuation = zero

Sample card for first-order hunting (multi-parameter continuation):

::

    Continuation = hfirst

Sample card for any LOCA continuation routine:

::

    Continuation = loca

**Technical Discussion**

Continuation refers to solving a given problem at a series of steps in which one of more numerical values (other than time) appearing in the governing equations and/or boundary conditions is varied between specified limits at specified step increments. If any LOCA algorithm is to be used, refer to the next paragraph; otherwise, set as follows: When there is only one such value, or continuation parameter, option zero or first is chosen and the information required to identify the parameter and its values is provided in the following cards in this section. When there are two or more parameters which are to be simultaneously varied using hunting continuation, option hzero or hfirst is chosen and the relevant information for each parameter is provided in the Hunting Specifications section.

To use any of the continuation routines in LOCA, option loca is chosen and the inputs required in this section will depend on the chosen routine. LOCA offers a number of continuation methods which can also be used with multiple parameters independently of the hunting algorithms.

The zero-order algorithms use the converged solution from each step as the initial guess for the next step, while the first-order algorithms perform a resolve to obtain a parameter sensitivity vector and use this vector to estimate the solution at the next step.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
