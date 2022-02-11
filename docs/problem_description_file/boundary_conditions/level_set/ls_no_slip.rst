**************
**LS_NO_SLIP**
**************

::

	BC = LS_NO_SLIP PF <integer>

-----------------------
**Description / Usage**
-----------------------

**(EMB/VECTOR MOMENTUM)**

This boundary condition is used to enforce the fluid/solid kinematic constraint for 1)
overset grid applications in which 2) the solid material is assumed rigid and therefore
has no internal stresses. It requires 3) a slaved phase function be defined along with 4)
a vector field of Lagrange multipliers.

A description of the input parameters follows:

============== =============================================================
**LS_NO_SLIP** Name of the boundary condition.
**PF**         This string indicates that this boundary condition is going to
               be applied along the zero contour of an embedded phase
               function (PF) field.
<integer>      This integer identifies the specific phase function field that
               is defining the contour. At the present time this integer
               should always be one.
============== =============================================================

------------
**Examples**
------------

An example:
::

   BC = LS_NO_SLIP PF 1

-------------------------
**Technical Discussion**
-------------------------

This boundary condition is used in the context of Goma’s overset grid capability. A
thorough treatment of this method is provided in the Goma document (GT-026.3) and
the user is directed there. However, a brief discussion of the nature of this boundary
condition is in order at this point.

The overset grid capability is used in problems in which a solid material is passing
through a fluid material. The solid and fluid materials both have there own meshes. In
the general problem, stresses and velocities must be transferred between each phase
and therefore there is two-coupling of the respective momentum and continuity
equations. This boundary condition, however, is used in the restricted case in which
the solid material is assumed to be rigid and having a prescribed motion. Therefore, the
coupling only proceeds in one direction : solid to fluid.

This boundary condition concerns itself with enforcing the kinematic constraint:

.. figure:: /figures/204_goma_physics.png
	:align: center
	:width: 90%

between the solid material with prescribed motion, :math:`\underline {\dot{x}}` , and the fluid whose velocity is, :math:`\gamma`. This kinematic constraint represents a new set of equations in the model for which
unknowns must be associated. In this case, we introduce a Lagrange multiplier vector
field, :math:`\underline{\gamma}` , at each node in the mesh. For fluid elements that do not intersect the fluid/
solid interface, these Lagrange multipliers are identically zero. They are non zero only
for those fluid elements that are crossed by the fluid/solid boundary. These Lagrange
multiplier fields couple the influence of the solid material on the fluid through body
force terms in the fluid momentum equations of the form:

.. figure:: /figures/205_goma_physics.png
	:align: center
	:width: 90%

When applying this boundary condition it is necessary to include Lagrange multiplier
equations equal to the number of dimensions in the problem. These are specified in the
equation section of the input deck. The shape and weight functions for these fields are
generally simple P0 functions. If one were to vector plot the components of the
Lagrange multiplier components, you get a general picture of the force interaction field
between the liquid and solid. This is sometimes informative.

A slaved phase function field is used to imprint the contour of the solid material on the
liquid mesh. The zero contour of this function is then used to evaluate the above line
integral. This phase function field is slaved to the solid material and is not evolved in
the conventional sense. Nonetheless, a single phase function field equation must be
included with the set of equations solved. In the phase function parameters section of
the input deck, the user must indicate that this phase function is slaved and also must
identify the sideset number of the boundary on the solid material which is the fluid/
solid interface.



--------------
**References**
--------------

No References. 

..
	TODO - Lines 60 and 72 have pictures that needs to be exhanged with equations.
