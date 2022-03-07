**************
TABLE_WICS
**************

::

	BC = TABLE_WICS SS <bc_id> {abscissa} {ordinate} {scale} {interpolation} [FILE = <fname>]

-----------------------
Description / Usage
-----------------------

**(WIC/VARIED)**

This boundary allows the user to supply boundary data for scalar weak integrated
boundary conditions. See the *TABLE_WICV* card for vector weak integrated boundary
conditions. A prime example of the use of the *TABLE_WICS* card is application of heat
flux for a thermal problem.

Definitions of the input parameters are as follows:

TABLE_WICS
    Name of the boundary condition (<bc_name>).                              
SS
    Type of boundary condition (<bc_type>), where **SS**
    denotes side set in the EXODUS II database.                              
<bc_id>
    The boundary flag identifier, an integer associated with
    *TABLE_WICS* that identifies the boundary location (side
    set in EXODUS II) in the problem domain.                                 
{abscissa}
    For one-dimensional tables (i.e. for use in 2D
    problems), the choices are restricted to one of the three
    coordinate directions. Use the strings X, Y or Z to
    identify the direction of choice. For two-dimensional
    tables (i.e. for use in 3D problems) use XY, XZ, YX,
    YZ, ZX, or ZY to denote the coordinate of the first and
    second columns in the table.                                             
{ordinate}
    This string identifies the equation of the weak integrated
    boundary term that the boundary data is added to. For
    example, use of the VELOCITY1 string will cause the
    table data to be used for the x-component of the liquid
    traction in the boundary integral for the x-momentum
    equation. See the following table.                                       

.. tabularcolumns:: |L|L|L|

+-------------------------------------------+--------------------------+------------------------------+
|**String**                                 |  **replaces**            |  **Equation**                |
+-------------------------------------------+--------------------------+------------------------------+
|VELOCITY1 or U                             |  liquid x-traction       |  R_MOMENTUM1                 |
+-------------------------------------------+--------------------------+------------------------------+
|VELOCITY2 or V                             |  liquid y-traction       |  R_MOMENTUM2                 |
+-------------------------------------------+--------------------------+------------------------------+
|VELOCITY3 or W                             |  liquid z-traction       |  R_MOMENTUM3                 |
+-------------------------------------------+--------------------------+------------------------------+
|TEMPERATURE                                |  diffusive energy flux   |  R_ENERGY                    |
+-------------------------------------------+--------------------------+------------------------------+
|MESH_DISPLACEMENT1 or DX or MESH_POSITION1 |  mesh x-traction         |  R_MESH1                     |
+-------------------------------------------+--------------------------+------------------------------+
|MESH_DISPLACEMENT2 or DY or MESH_POSITION2 |  mesh y-traction         |  R_MESH2                     |
+-------------------------------------------+--------------------------+------------------------------+
|MESH_DISPLACEMENT3 or DZ or MESH_POSITION3 |  mesh z-traction         |  R_MESH3                     |
+-------------------------------------------+--------------------------+------------------------------+
|SOLID_DISPLACEMENT1                        |  solid x-traction        |  R_SOLID1                    |
+-------------------------------------------+--------------------------+------------------------------+
|SOLID_DISPLACEMENT2                        |  solid y-traction        |  R_SOLID2                    |
+-------------------------------------------+--------------------------+------------------------------+
|SOLID_DISPLACEMENT3                        |  solid z-traction        |  R_SOLID3                    |
+-------------------------------------------+--------------------------+------------------------------+
|S[1-3][1-3]_[1-7]                          |  polymer mode traction   |  R_STRESS[1-3][1-3]_[1-7]    |
+-------------------------------------------+--------------------------+------------------------------+

{scale}
    A floating point scale multiplier which can be used to
    scale the tabular data. The boundary data used will be
    the product of {scale} and the tabular data.                             
{interpolation}
    This is the method chosen to interpolate between
    supplied data points.                                                    
                                                                                             
    For one-dimensional tables, the choices are **LINEAR**,
    which denotes linear interpolation, **QUADRATIC**,
    which denotes quadratic Lagrangian interpolation and
    requires an odd number of data points, and **QUAD_GP**,
    which denotes quadratic interpolation where the data
    points represent Gauss point values. 3N data points (see
    Technical Discussion) are required for **QUAD_GP** interpolation.        
                                                                                             
    For two-dimensional tables, **BIQUADRATIC** is
    currently the only choice. The first two columns of the
    table should define a rectangular, mapped grid where the
    second coordinate changes more quickly than the first.
    More complicated methods could be added latter.                          
[FILE = <fname>]
    The keyword "**FILE** =" indicates that the table data be
    read from a separate file identified by <fname>. This
    parameter is optional and if it is left out the table data
    will be read from the input deck itself following the
    *TABLE_WICS* card. In this latter case, the end of the
    table is signaled by the keywords "END TABLE". Note
    that the file specified by FILE = is fully *apreproable*,
    i.e., it will be preprocessed by APREPRO before
    reading if APREPRO is enabled.                                           

------------
Examples
------------

Following is a sample card:
::

     BC = TABLE_WICS SS 12 X TEMPERATURE QUADRATIC FILE =heatflux.table

:: 

    heatflux.table:
    
    0.0    1.0
    0.5    1.5
    1.0    1.75
    1.5    2.0
    2.0    2.0

-------------------------
Technical Discussion
-------------------------

The table data itself appears as columns of numbers. One-dimensional *TABLE_WICS*
tables have two columns (column1=abscissa, column2=ordinate), whereas twodimensional
*TABLE_WICS* tables have three columns (column1=abscissa1,
column2=abscissa2, column3=ordinate). *Goma* will try to read float values from any
line whose first parameter can be converted to a float.

The QUAD_GP interpolation option is meant for the case when the table data comes
from another finite element model or another *Goma* run and the data is most readily
available at the integration points of the finite element mesh. Hence, with quadratic
Gaussian quadrature, there are three data points per element. N is the number of
elements from the model that the data is coming from and therefore 3N data points are
the total expected.

The user is also referred to the section on **Boundary Condition Types** at the beginning
of the *Boundary Condition Specifications*. In particular, look at the discussion of
Weakly Integrated Conditions (WIC).

