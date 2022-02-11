*******
**MAT**
*******

::

	MAT = <char_string> <integer_list>

-----------------------
**Description / Usage**
-----------------------

This card represents the start of each material section in the *Problem Description File*. Thus, one MAT card is required for each material section. Definition of the input parameters are as follows:

=============== =========================================================
<char_string>   Filename of the material file from which all material
                properties for the current material will be read. The 
                material file’s name plus extension is char_string.mat, 
                and if the file is not present in the current working directory, the code will exit with the error message 
                “Not all Material Files found in current directory.”
<integer_list>  This is a list of space delimited integers that define the
                set of element blocks for which this material file is applicable; the integers are the element block ids 
                defined when the domain was meshed.
=============== =========================================================

------------
**Examples**
------------

The following specifies material file “sample.mat” applies to element blocks 1, 2, 3, 7,
and 9:
::

   MAT = sample 1 2 3 7 9

Note, the “.mat” extension is not specified explicitly, but appended to the character string by the code.

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

No References.