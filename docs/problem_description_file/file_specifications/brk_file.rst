*************
Brk File
*************

::

	Brk file = <file_name>

-----------------------
Description / Usage
-----------------------

<file_name>    
    A Brk file in the brk file syntax with specifications for material blocks

This optional card specifies the name of the Brk file for this problem, if one
does not exist goma will attempt to create one. The Brk file is used by the brk
utility to break the Exodus II files on parallel runs for each processor.

.. warning::
   Brk files can only be created on single processor runs.

------------
Examples
------------

Following is a sample card:
::

    Brk file = in.brk
