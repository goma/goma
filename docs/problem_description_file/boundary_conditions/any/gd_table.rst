************
GD_TABLE
************

::

	BC = GD_TABLE SS <bc_id> <equation_name> <integer1> <variable_name> <integer2> <scale> <interpolation> [FILE = <fname>]

-----------------------
Description / Usage
-----------------------

**(PCC/VARIED)**

This card is used to specify arbitrary, univariate (one abscissa and one ordinate: x1 - x2)
data for boundary conditions on two-dimensional boundaries, e.g., the inlet velocity
profile of a non-Newtonian fluid in a two-dimensional channel. The *GD_TABLE*
specification differs slightly from the other cards in this category: the data are scalable
and the data can be read from a file. Like the other *GD_** cards, this card can be used as
an additive building block for more complicated conditions. The examples below and
at the end of the *GD_** section will provide more detailed guidance.

Definitions of the input parameters are described next. Differences between this card
and other *GD_** cards are pointed out.

GD_TIME
    Name of the boundary condition (<bc_name>).
SS
    Type of boundary condition (<bc_type>), where **SS** denotes side set in
    the EXODUS II database.
<bc_id>
    The boundary flag identifier, an integer associated with <bc_type> that
    identifies the boundary location (side set in EXODUS II) in the problem
    domain.
<equation_name>
    A character string indicating the equation to which this boundary condition
    is applied. See the list of permissible values in the discussion above for
    *Category 1*. In contrast to other *GD_** cards, this parameter also serves
    to identify the equation that is being supplanted.
<integer1>
    Species number of the mass transport equation. The value should be 0 unless
    the <equation_name> is of type *R_MASS*.
<variable_name>
    A character string indicating the variable that should be used in the
    function. See the list of permissible values in the discussion above for
    *Category 1*. For this card, in contrast to other *GD_** cards, this
    parameter also identifies what value is to serve as abscissa when
    interpolating the table.
<integer2>
    Species number of the concentration variable.The value should be 0 unless
    the <variable_name> is of type *MASS_FRACTION*.
<scale>
    A floating point value by which to multiply the ordinate list after
    interpolation. It can be used to scale the table values or change their
    sign, *e.g*. C\ :sub:`0`, scale factor in f(x\ :sub:`1`) = C\ :sub:`0` * x\
    :sub:`2`
<interpolation>
    Specifies the method to use in interpolating between supplied data points.
    Currently the only choice available is *LINEAR*, which invokes a simple
    linear interpolation method. Alternative methods will/can be added latter
    as required or requested.

The table data will be read from within the input deck itself (following the GD_TABLE
BC card). The end of the table is signaled by the keywords "END TABLE." (See the
second example below.) An alternative to this method is to read a file with table data.

[FILE = <fname>]
    The optional keyword ‘*FILE* =’ indicates that the table data is to be read
    from a separate file identified by <fname>.

Note that this boundary condition card functions as every other GD condition, be it
*LINEAR, QUADRATIC, POLYNOMIAL*, or in this case *TABULAR*. It is used simple as
a piece of a residual on the appropriate equation. Hence, it usually requires more than
one GD card to completely specify the boundary condition.

------------
Examples
------------

Following is a sample card set in which the table data is to be read from an external file
called upstream_land.dat:
::

	BC = GD_LINEAR SS 1 R_MESH_NORMAL 0 MESH_POSITION2 0 0. -1.
	BC = GD_TABLE SS 1 R_MESH_NORMAL 0 MESH_POSITION1 0 1.0 LINEAR
	FILE=upstream_land.dat

This card set first creates a linear term in *MESH_POSITION2*, which is the y-coordinate
of the mesh points along side set 1. The second, *GD_TABLE* card then creates a table of
y-coordinate values based on x-mesh position. This boundary condition describes a
land/filet composite geometry with x-y data points.

Following is a sample card, where the table data is to be read directly from the input
deck:

::

        BC = GD_TABLE SS 1 R_MOMENTUM1 0 MESH_POSITION2 0 1.0 LINEAR

        $ r/R         Uz
        0.000000      1.666667
        0.050000      1.666458
        0.100000      1.665000
        0.150000      1.661042
        0.200000      1.653333
        0.250000      1.640625
        0.300000      1.621667
        .
        .
        0.900000      0.451667
        0.950000      0.237708
        1.000000      0.000000
        END TABLE

This table is used to specify the radial dependence of an axial velocity profile along the
specified side set.

-------------------------
Technical Discussion
-------------------------

This capability is widely used for geometry and velocity profile boundary conditions
that do not have a convenient closed form. Note that for geometry specifications you
cannot specify multi-valued functions, like for a cutback angle.



--------------
References
--------------

GTM-021.0: Multiparameter continuation and linear stability analysis on highly
deformable meshes in Goma, M. M. Hopkins, June 22, 2000


