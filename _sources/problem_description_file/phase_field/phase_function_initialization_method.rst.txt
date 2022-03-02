****************************************
**Phase Function Initialization Method**
****************************************

::

	Phase Function Initialization Method = {method_name} {parameter list}

-----------------------
**Description / Usage**
-----------------------

This card specifies the means by which the phase functions are initialized. After the
initial instance, subsequent instances of {model_name} {parameter_list} are used to
describe initializations of phase fields 2 through 5. This card constructs from a
representation of the starting interface shape, a value for the distance function at every
node in the mesh. The syntax of the card is as follows:

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

==========================================  =====================================================
**Projection**                              This method computes the initial phase function  
                                            field by
                                            calling a user-specified routine which returns the 
                                            signed
                                            distance function for a given point. It has no 
                                            parameter
                                            list after its name.
**Exodus**                                  Using this card indicates that the initial phase 
                                            function
                                            field is to be read from the exodus file specified 
                                            earlier
                                            (see *FEM file* and *Initial Guess* cards for **
                                            read_exoII**
                                            option). This card has no parameter list after its 
                                            name.
**Nodeset** <integer1> **EB** <integer2>    This method establishes the initial location of the
                                            interface as the boundary between two element blocks.
                                            The value <integer1> is the nodeset identification
                                            number for an internal nodeset defined to exist at 
                                            the
                                            interface between the two element blocks. The 
                                            character
                                            string **EB** is required. The integer <integer2> is 
                                            the
                                            element block id number to which positive values of
                                            phase function function is going to be assigned.
**Surfaces** <integer>                      This card establishes the initial phase function   
                                            function
                                            by referring to a set of primitive geometric 
                                            objects. It is
                                            the easiest to use and the most general. The integer
                                            value <integer> is the number of surface objects 
                                            that are
                                            used to construct the initial interface. This number 
                                            of
                                            **SURF** object cards must follow this card. This is 
                                            the
                                            syntax of the **SURF** object card:

                                            **SURF** = {object_name} {float list}  

                                            {object_name}: a character string identifying the
                                            type of geometric object. Options are: **PLANE,**
                                            **CIRCLE, SPHERE,** SS, **USER.**

                                            {float list}: geometric parameters associated with
                                            each object as float values.
==========================================  =====================================================

The following is the syntax and description for each geometric
object option, i.e., the “{object_name} {float list}” part of **SURF**

.. tabularcolumns:: |l|L|

==========================================  ====================================================
**PLANE** <nx. <ny> <nz> <d>                This card constructs a planar interface surface.The 
                                            float
                                            values <nx>, <ny>, <nz> define a vector normal to 
                                            this
                                            plane with the restriction that the sign of the 
                                            vector must
                                            be such that it points from the negative side of the
                                            interface to the positive side of the interface. 
                                            The 
                                            float value <d> effectively represents the distance 
                                            of the
                                            plane from the origin. Its value must be set, 
                                            however, so
                                            that the dot product of any position vector to a 
                                            point on
                                            the desired plane and the vector (nx,ny,nz) must be
                                            equal to <d> (it is a property of planes that this 
                                            number
                                            is independent of the point on the plane that is 
                                            chosen).
**CIRCLE** <cx> <cy> <radius>               This card constructs a circular interface surface 
                                            in a
                                            two-dimensional domain. The float values <cx> <cy>
                                            identify the coordinates of the center of the 
                                            circle. The
                                            float value <radius> establishes the radius of the 
                                            curve.
                                            By definition, points interior to the circle are 
                                            assigned
                                            negative phase function function values.
**SPHERE** <cx> <cy> <cz> <radius>          This card constructs a spherical interface surface 
                                            in a
                                            three-dimensional domain. The float values <cx> <cy>
                                            <cz> identify the coordinates of the center of the 
                                            circle.
                                            The float value <radius> establishes the radius of 
                                            the
                                            sphere. By definition, points interior to the 
                                            sphere are
                                            assigned negative phase function function values.
SS {ss_id}                                  This card uses an existing sideset in the problem 
                                            as a
                                            defined geometric object for construction of an
                                            interface. The parameter <ss_id> identifies this 
                                            sideset.
**USER** {user-defined float list}          This card indicates the user has defined an object
                                            function using the supplied parameter float list 
                                            that
                                            returns a signed distance value when supplied with 
                                            the
                                            coordinates of a point in space. This object 
                                            function
                                            should appear in the function call *user_init_object
                                            * in the file **user_pre.c.**
**SM_object** {object_type} {object_name}   This card allows the user to initialize the phase 
                                            function
                                            location by using a piece of solid model geometry. 
                                            The
                                            solid model object_type can be either **FACE** or **
                                            BODY.**
                                            A 2D initialization uses the boundary of the 
                                            specified
                                            FACE (or surface) as the 0 phase function. A 3D 
                                            initialization uses the boundary of the specified 
                                            BODY (or volume) as the 0 phase function.
==========================================  ====================================================

------------
**Examples**
------------

Three examples of initialization methods for a single phase function are provide below:
::

	Phase Function Initialization Method = Nodeset 20 EB 1

::

	Phase Function Initialization Method = Surfaces 3
		SURF = PLANE -1. 0. 0. -3.
		SURF = CIRCLE -2 0 1
		SURF = CIRCLE -3 0 0.5

::

	Phase Function Initialization Method = SM_object BODY my_blob

-------------------------
**Technical Discussion**
-------------------------

Please consult Level Set Initialization Method card for discussion.

--------------
**References**
--------------

GT-020.1: Tutorial on Level Set Interface Tracking in GOMA, February 27, 2001, T.A.
Baer
