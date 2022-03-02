***********************************
Level Set Initialization Method
***********************************

::

	Level Set Initialization Method = {method_name} {parameter list}

-----------------------
Description / Usage
-----------------------

This card specifies the means by which the level set function is initialized. That is, it
constructs from a representation of the starting interface shape, a value for the distance
function at every node in the mesh. The syntax of the card is as follows:

{method_name}
    A character string which identifies the initialization option desired.
    Choices for this string are: **Projection, Exodus, Nodeset, Surfaces,
    SM_object.**

{parameter list}
    This is a variable parameter list specific to each option. The nature of it
    for each method is detailed in the syntax descriptions below.

Below are the exact syntax used for each initialization method, a brief description of
the method and a specification of any additional required parameters.

.. tabularcolumns:: |l|L|
 
========================================  ============================================================
**Projection**                             This method computes the initial level set field by
                                           calling a user-specified routine which returns the signed
                                           distance function for a given point. It has no parameter
                                           list after its name.
**Exodus**                                 Using this card indicates that the initial level set 
                                           field is
                                           to be read from the exodus file specified earlier (see
                                           *FEM file* and *Initial Guess* cards for **read_exoII**
                                           option). This card has no parameter list after its name.
**Nodeset** <integer1> **EB** <integer2>   This method establishes the initial location of the
                                           interface as the boundary between two element blocks.
                                           The value <integer1> is the nodeset identification
                                           number for an internal nodeset defined to exist at the
                                           interface between the two element blocks. The character
                                           string **EB** is required. The integer <integer2> is the
                                           element block id number to which positive values of
                                           level set function is going to be assigned.
**Surfaces** <integer>                     This card establishes the initial level set function by
                                           referring to a set of primitive geometric objects. It is the
                                           easiest to use and the most general. The integer value
                                           <integer> is the number of surface objects that are used
                                           to construct the initial interface. This number of **SURF**
                                           object cards must follow this card. This is the syntax of
                                           the **SURF** object card:

                                           SURF = {object_name} {float list}

                                           {object_name}: a character string identifying the
                                           type of geometric object. Options are: **PLANE**,
                                           **CIRCLE,** **SPHERE,** SS, **USER.**

                                           {float list}: geometric parameters associated with
                                           each object as float values
========================================  ============================================================

The following is the syntax and description for each geometric
object option, i.e., the “{object_name} {float list}” part of **SURF**

.. tabularcolumns:: |l|L|

========================================= ============================================================
**PLANE** <nx. <ny> <nz> <d>               This card constructs a planar interface surface. The float
                                           values <nx>, <ny>, <nz> define a vector normal to this
                                           plane with the restriction that the sign of the vector must
                                           be such that it points from the negative side of the
                                           interface to the positive side of the interface. The float
                                           value <d> effectively represents the distance of the
                                           plane from the origin. Its value must be set, however, so
                                           that the dot product of any position vector to a point on
                                           the desired plane and the vector (nx,ny,nz) must be equal to 
                                           <d> (it is a property of planes that this number
                                           is independent of the point on the plane that is chosen).
**CIRCLE** <cx> <cy> <radius>              This card constructs a circular interface surface in a
                                           two-dimensional domain. The float values <cx> <cy>
                                           identify the coordinates of the center of the circle. The
                                           float value <radius> establishes the radius of the curve.
                                           By definition, points interior to the circle are assigned
                                           negative level set function values.
**SPHERE** <cx> <cy> <cz> <radius>         This card constructs a spherical interface surface in a
                                           three-dimensional domain. The float values <cx> <cy>
                                           *<cz>* identify the coordinates of the center of the circle.
                                           The float value <radius> establishes the radius of the
                                           sphere. By definition, points interior to the sphere are
                                           assigned negative level set function values.
SS {ss_id}                                 This card uses an existing sideset in the problem as a
                                           defined geometric object for construction of an
                                           interface. The parameter <ss_id> identifies this sideset.
**USER** {user-defined float list}         This card indicates the user has defined an object
                                           function using the supplied parameter float list that
                                           returns a signed distance value when supplied with the
                                           coordinates of a point in space. This object function
                                           should appear in the function call *user_init_object* in the
                                           file **user_pre.c.**
**SM_object** {object_type} {object_name}  This card allows the user to initialize the level set
                                           location by using a piece of solid model geometry. The
                                           solid model object_type can be either **FACE** or **BODY.**
                                           A 2D initialization uses the boundary of the specified
                                           FACE (or surface) as the 0 level set. A 3D initialization
                                           uses the boundary of the specified BODY (or volume)
                                           as the 0 level set.
========================================= ============================================================

------------
Examples
------------

Two examples of initialization methods are provide below:
::

	Level Set Initialization Method = Nodeset 20 EB 1

::

	Level Set Initialization Method = Surfaces 3
            SURF = PLANE -1. 0. 0. -3.
		SURF = CIRCLE -2 0 1
		SURF = CIRCLE -3 0 0.5

::

	Level Set Initialization Method = SM_object BODY my_blob

-------------------------
Technical Discussion
-------------------------

The **Projection** initialization method was developed early in the level set
development process. It has since been superseded by other more easily used
methods. It is still supported primarily for the use of developers. Users wanting a
complicated interface shape for which they can supply an appropriate distance
function should user the USER surface object option under the Surfaces
initialization method.

The **Exodus** method deserves little comment. It should be used when restarting
level set computations from a preexisting solution.

The **Nodeset** method allows the user to make use of the sophisticated solid body
manipulation software in meshing packages like CUBIT. The procedure for using
this method is to create a domain which contains two element blocks. The desired
starting point for the interface should lie on the curve or surface which these two
blocks have in common. A single nodeset should be defined over this entire curve
or surface. The nodeset identification number should be the first integer parameter
specified on the card. Also note that one of the blocks must be designated as the
“positive” block. This means then when initialized the values of the level set
function in this block will be positive. The values in the other block will be
negative. Note that this initialization method can only by used for problems that
have exactly two blocks, no more.

The **Surfaces** initialization method is the most useful method for initialization. It
draws from the fact that it is relatively easy to determine the distance to simple
geometric objects (planes, circles, spheres, etc.). Further, it permits initialization
using more than one of these objects so that relatively complicated initial interface
locations can be constructed. However, the user should recognize that this method
is still somewhat unsophisticated in its approach so there are some caveats
associated with its use. The primary point is that surface objects should never
intersect anywhere within the domain of interest, otherwise it is more than likely
that the starting interface shape will not be what the user expects.

The **SM_object** initialization method allows the user to use solid model geometry
to initialize 2D and 3D level sets. Certain 2D geometries can be created using only
Goma input commands (see *FACE*). Other 2D geometries, and all 3D geometries,
can be accessed via an ACIS .sat file. The usual way to do this is for the user to
create their desired geometry within Cubit (or, import solid model geometry from
elsewhere into Cubit). Faces (or surfaces) should be created for 2D initialization,
and bodies (or volumes) should be created for 3D initialization. The *boundary* of
the object is used to initialize the level set. The geometry should be named within
Cubit and exported to an ACIS .sat file via Cubit’s export acis
“filename” ascii command. This same file should be read in via the *ACIS
file* command in the Geometry Specifications section. The solid model geometry is
then available for the *Level Set Initialization Method* command. (Note that the
Geometry Specifications section usually comes after the *Level Set Initialization
Method* command; this is OK).

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
