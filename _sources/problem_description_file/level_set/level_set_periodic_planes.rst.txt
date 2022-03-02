*****************************
Level Set Periodic Planes
*****************************

::

	Level Set Periodic Planes = <float1> <float2> <float3> <float4> <float5> <float6>

-----------------------
Description / Usage
-----------------------

This card directs the level-set renormalization to accommodate peroidic boundary
conditions (see Advanced Capabilities Manual AC_Periodic capability). The periodic
boundary conditions on the level set field are not compatible with renormalization
unless this capability is specified.

<float1>
    x-coordinate value of first periodic boundary.

<float2>
    x-coordinate value of second periodic boundary. If equivalent to float1
    than this direction is not periodic.

<float3>
    y-coordinate value of first periodic boundary.

<float4>
    y-coordinate value of second periodic boundary. If equivalent to float3
    than this direction is not periodic.

<float5>
    y-coordinate value of first periodic boundary.

<float6>
    y-coordinate value of second periodic boundary. If equivalent to float5
    than this direction is not periodic.

------------
Examples
------------

Two examples of initialization methods are provide below:
::

	Level Set Periodic Boundary = -0.5 0.5 0 0 0 0

This card instructs renormalization to accommodate the x-directed-boundaries to be
considered as periodic relative to the level-set field.

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
