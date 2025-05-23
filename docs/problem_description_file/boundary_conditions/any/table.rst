*********
TABLE
*********

::

	BC = TABLE SS <bc_id> {X|Y|Z|TIME} {ordinate} [species] {interpolation} [FILE = <fname>] [NAME = <identifier>]

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED))**

This boundary condition is a stand-alone version of the more complicated *GD_TABLE*
card. It allows the user to supply arbitrary univariate (one abscissa and one ordinate)
data about the spatial variation of unknowns fields on a boundary. The abscissa will be
one of the three spatial coordinates or time and the ordinate is one of a choice of
unknown field variables. All *TABLE_BC* conditions must have attached tabular data as
a list of paired float values either directly following the card or in a separate file
(identified on the card). The list of data pairs is terminated by the string “END TABLE”
on its own line.

Definitions of the input parameters are as follows:

TABLE
    Name of the boundary condition.
SS
    Type of boundary condition (<bc_type>), where **SS** denotes side set in
    the EXODUS II database.
<bc_id>
    The boundary flag identifier, an integer associated with <bc_type> that
    identifies the boundary location (side set in EXODUS II) in the problem
    domain.
{X|Y|Z|TIME}
    A char_string that identifies the independent table variable (abscissa).
    The strings X,Y, and Z refer of course to the three spatial coordinates.
    Depending on the choice here, the x, y, or z coordinate value at a given
    point, respectively, is used to obtain an interpolated ordinate value using
    the attached table data. If the TIME string appears here, however, the
    current simulation time is used to interpolate an ordinate value. This
    single value is applied uniformly to the sideset.
{ordinate}
    This string associates a variable type with the values of the ordinate in
    the attached table. It also identifies the equation that is supplanted by
    the boundary condition on the sideset. The following table lists the
    available string choices and the corresponding equation component clobbered
    by the boundary condition.

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
VELOCITY1 or U                             R_MOMENTUM1
VELOCITY2 or V                             R_MOMENTUM2
VELOCITY3 or W                             R_MOMENTUM3
MASS_FRACTION or Y or SPECIES              R_MASS
TEMPERATURE                                R_ENERGY
MESH_DISPLACEMENT1 or DX                   R_MESH1
MESH_DISPLACEMENT2 or DY                   R_MESH2
MESH_DISPLACEMENT3 or DZ                   R_MESH3
PRESSURE or P                              R_PRESSURE
SOLID_DISPLACEMENT1 or DX_RS               R_SOLID1
SOLID_DISPLACEMENT2 or DY_RS               R_SOLID2
SOLID_DISPLACEMENT3 or DZ_RS               R_SOLID3
SHEAR_RATE or SH                           R_SHEAR_RATE
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11                                        R_STRESS11
S12                                        R_STRESS11
S22                                        R_STRESS11
S13                                        R_STRESS11
S23                                        R_STRESS11
S33                                        R_STRESS11
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_1                                      R_STRESS11_1
S12_1                                      R_STRESS12_1
S22_1                                      R_STRESS22_1
S13_1                                      R_STRESS13_1
S23_1                                      R_STRESS23_1
S33_1                                      R_STRESS33_1
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_2                                      R_STRESS11_2
S12_2                                      R_STRESS12_2
S22_2                                      R_STRESS22_2
S13_2                                      R_STRESS13_2
S23_2                                      R_STRESS23_2
S33_2                                      R_STRESS33_2
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_3                                      R_STRESS11_3
S12_3                                      R_STRESS12_3
S22_3                                      R_STRESS22_3
S13_3                                      R_STRESS13_3
S23_3                                      R_STRESS23_3
S33_3                                      R_STRESS33_3
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_4                                      R_STRESS11_4
S12_4                                      R_STRESS12_4
S22_4                                      R_STRESS22_4
S13_4                                      R_STRESS13_4
S23_4                                      R_STRESS23_4
S33_4                                      R_STRESS33_4
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_5                                      R_STRESS11_5
S12_5                                      R_STRESS12_5
S22_5                                      R_STRESS22_5
S13_5                                      R_STRESS13_5
S23_5                                      R_STRESS23_5
S33_5                                      R_STRESS33_5
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_6                                      R_STRESS11_6
S12_6                                      R_STRESS12_6
S22_6                                      R_STRESS22_6
S13_6                                      R_STRESS13_6
S23_6                                      R_STRESS23_6
S33_6                                      R_STRESS33_6
============================= ============ ================================

.. tabularcolumns:: |L|l|L|

============================= ============ ================================
**String**                    **replaces** **Equation***
----------------------------- ------------ --------------------------------
S11_7                                      R_STRESS11_7
S12_7                                      R_STRESS12_7
S22_7                                      R_STRESS22_7
S13_7                                      R_STRESS13_7
S23_7                                      R_STRESS23_7
S33_7                                      R_STRESS33_7
============================= ============ ================================

[species]
    An optional integer parameter that identifies the index of the appropriate
    species. Note, it should appear only when the <ordinate> string is
    *MASS_FRACTION*.
{interpolation}
    A char_string parameter that identifies the method chosen to interpolate
    between the attached table data points. For one-dimensional tables, the
    choices are *LINEAR*, which denotes simple linear interpolation, and
    *QUADRATIC*, which denotes quadratic Lagrangian interpolation. Note that
    the latter requires an odd number of data points be supplied in the table.
[FILE = <fname>]
    The optional char_string keyword "**FILE** =" indicates that the table data
    be read from a separate file identified by <fname>. This parameter is
    optional and if it is left out the table data will be read from the input
    deck itself following the *TABLE BC* card. Note that the file specified by
    <fname> will be first preprocessed by APREPRO if that option was enabled on
    the command line. This is a useful feature that allows for a quick way to
    introduce analytic expressions onto boundaries.
[NAME = <identifier>]
    The optional char_string keyword *NAME* = allows for a set of table data to
    be attached to the char_string parameter <identifier>. This option can only
    be used if the table data is read from a separate file identified by *FILE*
    = <*fname*>. In this case, the file <fname> is scanned for the char_string
    “identifier:” (note the colon). Once found the table data is read until
    encountering *END TABLE*. This option permits multiple sets of data in the
    same file.

The second half of the *TABLE_BC* implementation is the tabular data itself. In the
*TABLE* boundary condition, it consists of a set of paired float values, each pair on its
own line. This data should follow directly after the *TABLE* boundary condition card if
the *FILE* = option is not used. If a value for <fname> is supplied, the table data should
be written in the file so indicated. Note that in most implementations of UNIX,
<fname> can include a complete path specification in case the datafile is in a different
directory than the run directory. In either case, input deck or separate file, the set of
data table pairs should always be terminated by the string *END TABLE* to terminate
reading of the data. When reading the table data, *Goma* attempts to read a float value on
each line. If it is unsuccessful, e.g., a string might start the line, it will proceed to the
next line. If it is successful, it will attempt to read a second float value to complete the
data pair. An unsuccessful read here is an error. Once the second value is read,
however, the remainder of the line is discarded and the next line is read. This procedure
permits inclusion of comments within. See the next section for some examples.

Thus,

::

        3. 1.e-4
        1. 3. % this is a good example
        $ 1. 40.0
        $ I have no idea where the following data came from
            3.4   2.1
            1.e-2   6000.0

will result in four data points being read, whereas, both of the following

::

        6.443   3.43c
        5.4099   % 099.0

will result in an error.

------------
**Examples**
------------

The following is an example of a tabular data set that will be read correctly
::

        $ This data came from M. Hobbs. God only knows where he got it.
        T   k
        0.5 1.e-4
        1. 15.   % I’m not particularly sure about this one.
        3.4   8.1
        5.6   23.0
        $ 1.0 40.0

In this case, four data pairs will be read to form the table.

Example usage of the *TABLE* card follows:

        * Setting the u-velocity on an inlet boundary for a power law fluid:

::

        BC = TABLE SS 1 Y U LINEAR
        $ r/R Ux
        0.000000 1.666667
        0.050000 1.666458
        0.100000 1.665000
        0.150000 1.661042
        0.200000 1.653333
        0.250000 1.640625
        0.300000 1.621667
        ..
        ..
        0.900000 0.451667
        0.950000 0.237708
        1.000000 0.000000
        END TABLE

|

        * Setting the inlet concentration profiles for species 0 and species 1 from data in
          y.table:

::

        BC = TABLE SS 1 Y SPECIES 0 QUADRATIC FILE = y.table NAME = y0
        BC = TABLE SS 1 Y SPECIES 1 QUADRATIC FILE = y.table NAME = y1

|

        * The file y.table contains:

::

        y0:
                0.   1.0
                0.25 0.75
                0.5 0.60
                0.75 0.30
                1.0 0.20
        END TABLE
        y1:
                0. 0.0
                0.25 0.2
                0.5 0.3
                0.75 0.5
                1.0 0.8
        END TABLE

|

        * Setting a temperature history on a sideset

::

        BC = TABLE SS 1 TIME TEMPERATURE LINEAR
        0.0   0.0
        10.0   373.0
        40.0   373.0
        50.0   500.0
        100.0   500.0
        150    0.0
        100000.0   0.0
        END TABLE

-------------------------
**Technical Discussion**
-------------------------

The *TABLE* boundary condition provides similar functionality to the *GD_TABLE*
boundary condition but with a simplified interface the notion behind both cards is that
often information on boundaries is known only as a set of data points at specific
positions on the boundary. The *TABLE* boundary condition can use that boundary
information to provide interpolated values at nodal locations and then impose them as a
strong point collocated condition.

Interpolation orders for this method are limited to *LINEAR* and *QUADRATIC* with the
latter requiring an odd number of data points be supplied in the table.
