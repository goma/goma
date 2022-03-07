************
FEM File
************
::

	FEM file = <file_name>

-----------------------
Description / Usage
-----------------------

This required card specifies the name of the EXODUS II finite element mesh file. Any
EXODUS II file name is permissible, as specified below.

.. tabularcolumns:: |l|L|

<file_name>          
    A file name of the form *prefix*.exoII. The *prefix* portion is any
    user-specified alpha-numeric string, which can be used as a problem-type
    descriptor. Preprocessors and postprocessors (like AVS) might require the
    “.exoII” suffix so it is a required part of the file designation. The
    maximum length of the file name is 85 characters.

------------
Examples
------------

Following is a sample card:
::

	FEM file = in.exoII

-------------------------
Technical Discussion
-------------------------

This file contains the finite element discretization of the problem domain. Finite
element mesh files from other preprocessors may be used with *Goma* as long as a
translator from the preprocessor’s output format to the EXODUS II format is available
to the analyst.

--------------
References
--------------

The EXODUS II format is documented in:

* EXODUS II: A Finite Element Data Model, Schoof, L. A. and V. R. Yarberry, SAND92-2137, Sandia National Laboratories, Albuquerque, NM.
