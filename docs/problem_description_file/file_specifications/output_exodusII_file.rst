*************************
Output EXODUS II File
*************************

::

	Output EXODUS II file = <file_name>

-----------------------
Description / Usage
-----------------------

This required card specifies the name of the output EXODUS II file. Any EXODUS II
file name is permissible, as specified below.

<file_name> 	
    A file name of the form ``*prefix*.exoII``. The *prefix* portion is any
    user-specified alpha-numeric string, which can be used as an output file
    descriptor.

This EXODUS II file contains a replica of the input mesh and boundary condition
information exactly as it was provided in the *FEM file*, but has appended to it the
solution field information appropriate to the problem type. If the name of this output
EXODUS II file <file_name> is identical to the name of the input EXODUS II file (as
specified in the *FEM file* card), then no replication of the input mesh data is performed
and any results are simply appended to it.

------------
Examples
------------

Following is a sample card:
::

	Output EXODUS II file = out.exoII

-------------------------
Technical Discussion
-------------------------

Although allowed, it is not advisable to make this file name the same as the file name
input on the *FEM file* card.


--------------
References
--------------

The EXODUS II format is documented in:

* EXODUS II: A Finite Element Data Model, Schoof, L. A. and V. R. Yarberry, SAND92-2137, Sandia National Laboratories, Albuquerque, NM.
