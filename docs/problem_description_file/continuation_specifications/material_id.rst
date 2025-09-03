Material id
-----------------

**Material id** = <integer>

**Description/Usage**

This card is required for continuation problems when the specified Continuation Type is either "MT" or "UM". It identifies which of the materials specified for the problem has the property which will be used as the continuation parameter, and corresponds to the block number assigned to this material in the input Exodus file. If there is only one material, this number will always be 1. If there are two or more materials, they would be numbered starting with 1.

**Examples**

For a problem with only one material, use:

::

    Material id = 1

If there are several materials and the continuation parameter is a property of aluminum (with properties given in file "aluminum.mat") which occupies block 4 of the input Exodus file, there will be a subsection for aluminum in the Problem Description section which will start with the card:

::

    MAT = aluminum 4

In this case, use:

::

    Material id = 4

**Technical Discussion**

This index number differs from most others in that it is one-based (to be consistent with the input file) even though the Goma internal index is zero-based. Accordingly, there is never a case where material 0 would be specified.

**Theory**

No Theory.

**FAQs**

No FAQs.

**References**

No References.
