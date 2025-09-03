Boundary condition ID
----------------------------

**Boundary condition ID** = <integer>

**Description/Usage**

This card is required for continuation problems when the specified Continuation Type is BC or AC.

<integer>
    the zero-based numerical position of the BC or AC card (within the list in the Boundary or Augmenting Condition Specifications section) which contains the continuation parameter. If there are N BC or AC cards, the top card would be number 0 and the bottom card would be number N-1.

**Examples**

If the continuation parameter is included in the third BC card from the top of the list, this card should be:

::

    Boundary condition ID = 2

If the continuation parameter is included in the first (or only) AC card, use:

::

    Boundary condition ID = 0

**Technical Discussion**

Note that the BC list written to the screen when Goma is launched with the -bc_list command line option (which prints the active boundary conditions with their assigned indices) is one-based, such that the above BC would be numbered 3. If this list is used to determine the ID number, be careful to subtract one from the number on the screen.

Note that the BC ID cards for continuation are being "borrowed" for similar use by augmenting conditions, in lieu of creating a second set of cards.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References
