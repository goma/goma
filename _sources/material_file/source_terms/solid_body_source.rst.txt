*********************
**Solid Body Source**
*********************

::

   Solid Body Source = CONSTANT <species_number> <float1> <float2> <float3>

-----------------------
**Description / Usage**
-----------------------

This card is used to specify the model for the body force source term on the solid
mechanics momentum equations. This card is used most to impose gravitational forces
on solid phase material elements in the problem. It can be also used to impose body
forces on the pseudo-solid mesh material if that is desirable.

Definitions of the input parameters are as follows:

+--------------------------+-------------------------------------------------------------------------------------+
|**CONSTANT**              |A string identifying the constant force model. Currently, this is the only body force|
|                          |model for solid materials.                                                           |
|                          |                                                                                     |
|                          | * <float1> - the x-component of the body force vector in (F/L3).                    |
|                          | * <float2> - the y-component of the body force vector in (F/L3).                    |
|                          | * <float3> - the z-component of the body force vector in (F/L3).                    |
+--------------------------+-------------------------------------------------------------------------------------+
|**JXB**                   |A string identifying a body force model based on external current density fields J   |
|                          |and external magnetic fields B. See Technical discussion.                            |
|                          |                                                                                     |
|                          | * <float1> - a scale factor, usually set to 1.                                      |
+--------------------------+-------------------------------------------------------------------------------------+


------------
**Examples**
------------

The following is a sample input card:

::

   Solid Body Source = CONSTANT 0.0 0.0 -2000.0

-------------------------
**Technical Discussion**
-------------------------

Just as there is a body force vector that can be applied to fluid material regions, there is
a capability to apply a similar body force term to solid material regions. Most often this
is used to apply gravitational forces in which case the values of the components
supplied on this card would be the solid density multiplied by the gravitational
acceleration vector.

The **JXB** model requires external nodal fields loaded through External Field capability.
These fields must be named JE_N_1, JE_N_2, and JE_N_3 for the three components of
the current density and BE_N_1, BE_N_2, and BE_N_3 for the three components of
the magnetic field.



--------------
**References**
--------------

No References.