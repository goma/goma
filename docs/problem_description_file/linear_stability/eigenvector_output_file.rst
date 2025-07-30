Eigenvector output file
=======================

.. code-block:: none

    Eigenvector output file = <string>

Description/Usage
-----------------

This optional card is used to select a custom base name for the EXODUS II files which 
will be created for writing eigenvectors from the ARPACK eigensolver. The name 
should be entered as follows:

basename.exoII
    Filenames will be basename_modei.exoII

If this card is absent, the default name is LSA.exoII.

Examples
--------

Here is a sample card:

::

    Eigenvector output file = base-EV.exoII

Technical Discussion
--------------------

The actual filenames are formed from the base filename as follows:

For mode i (zero-based) when not using 3D of 2D LSA: basename_modei.exoII.

When using 3D of 2D LSA with wave number N: basename_modei_wn=N.exoII.

When requesting linear stability analysis during a LOCA continuation run, the 
eigenvectors from all continuation steps corresponding to a given mode number (or 
mode number/wave number combination for 3D of 2D) will be written sequentially to 
the same EXODUS II file. The time stamp for each step will be the eigenvalue, rather 
than the continuation parameter value. Complex eigenvectors are written as two steps 
in the same file: first the real part (time-stamped with the real part of the eigenvalue), 
then the imaginary part (time-stamped with the imaginary part of the eigenvalue).

The inputs Eigenvalue output frequency and Eigenvector output frequency are 
independent, so if a frequency m is specified for eigenvalues and n for eigenvectors, 
and n is not a multiple of m, then eigenvectors will be written on step numbers which 
are common multiples of m and n.

This card is not applicable to the eggroll eigensolver, which has a different default file 
naming convention with no options.

Theory
------

No Theory.

FAQs
----

No FAQs.

References
----------

No References.
