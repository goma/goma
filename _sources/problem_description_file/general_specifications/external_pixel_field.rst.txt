************************
**External Pixel Field**
************************

::

	External Pixel Field = <char_string1> {Q1|Q2} <file_name> <integer1>

-----------------------
**Description / Usage**
-----------------------

This optional card format provides a mechanism for reading-in pixel fields which are
converted (mapped, with a least squares algorithm) to finite element fields with the
chosen interpolation. After GOMA execution these fields are output in exodusII
format in a file map.exoII. Please see discussion below and Tutorial GT-038 for more
details and important tips.

<char_string1>     
    Name of the nodal field to be read; it should correspond to a nodal
    variable name you wish to have in the output EXODUS II file. If you
    subsequently wish to read the field in again and use as an EXTERNAL_FIELD
    model on other material property card, the chosen name matters.

{Q1 | Q2}          
    The type of interpolation to be applied to the external pixel field
    variable field. Possible values are as follows:

    * Q1  -  Linear
    * Q2  -  Quadratic

<file_name>        
    Name of the text file name with the pixel points. The pixel field format in
    this file should be as follows:

    * # pixel points
    * x_1, y_1, z_1 value
    * x_2, y_2, z_2 value
    * ...
    * x_N, y_n, z_n value

<integer1>         
    Material block ID to which the pixel field is mapped.

------------
**Examples**
------------

::

	External Pixel Field = HEIGHT Q1 tread.txt 1

-------------------------
**Technical Discussion**
-------------------------

Please consult the tutorial GT-038 before using this capability. Many user tips are
given together with a more thorough explanation on the proper use This capability is
extremely memory intensive, and excessive grid sizes and pixel densities can blow out
the memory on your machine. As of 12/22/2012 (the end of the Mayan calendar)
these fields are used typically to bring in pattern maps for scaling porous media and
lubrication height properties.

SAT, HEIGHT, PERM, CROSS_PERM, SH_SAT_CL_POROSITY, etc. These are
specially designated external fields which are mapped to variations in these properties
corresponding to thin porous media. Please see GT-038.


--------------
**References**
--------------

GT-038.0: Pixel-to-Mesh Tool Tutorial for GOMA. P R. Schunk, Memo to
distribution, 10 November 2009.

.. 
	TODO - Lines 35-29, need to be formatted in such a way that the correct message is being depicted. a table will no do this.
