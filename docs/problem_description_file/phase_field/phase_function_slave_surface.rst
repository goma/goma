********************************
**Phase Function Slave Surface**
********************************

::

	Phase Function Slave Surface = <char_string>

-----------------------
**Description / Usage**
-----------------------

This card is used to designate that the phase function degree of freedom is being slaved
to a boundary. This card is used primarily in the overset grid algorithm in which a
phase function field is slaved to the surface of the embedded body.

<char_string>
    YES|ON (not case sensitive) will allow the phase function field to be
    slaved to a surface. Currently, no support is given to more than one slaved
    function fields or to problems in which there are slaved and unslaved
    (free?) phase function fields.

------------
**Examples**
------------

A typical length scale input card looks like: 
::

	Phase Function Slave Surface = yes

-------------------------
**Technical Discussion**
-------------------------

One of the nice properties of level set/phase function fields is that they can be used to
find distances from surfaces. This function can be used quite apart from their abilities
to track interfaces. Including this card informs Goma that the phase function 1 field is
going to be used in this capacity and that no PDE is going to be solved to evolve it.
Instead, the values of this field will be “slaved” to a specific surface in the problem and
their values will be determined in reference to this surface in a process very reminicent
of renormalization.

The overset grid method makes use of a slaved phase function field. In that case, the
phase function field is slaved to the surface of the embedded object. As the embedded
object moves through the flow field, the slaved phase function values will be updated
by determining the distance of a given node to the object’s surface. This slaved phase
function field is then used in a variety of ways to compute the influence of the
embedded object on the flow and stresses of the surrounding fluid.

