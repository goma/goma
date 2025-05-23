************************************
Level Set Renormalization Method
************************************

::

	Level Set Renormalization Method = {char_string}

-----------------------
Description / Usage
-----------------------

This card indicates the method to be used to renormalize the level set function during
the course of the computation. The syntax of this card is as follows:

{char_string}
    A character string which specifies the type of method for renormalization.
    Choices for this string are: **Huygens, Huygens_Constrained, Correction.**

Each method is described below; see also the Technical Discussion.

Huygens
    In this method a set of m points P is constructed:

       :math:`\mathbf{P} = \left\{ \left( x_i, y_i, z_i \right), \,  
       i = 1,2, \ldots m | \quad \phi_j \left( x_i, y_i, z_i \right) 
       = 0 \right\}`

    which in a sense represent a discretization of the 
    interface location. The finite element interpolation
    functions are used to find exact locations for these
    points. For each mesh node :math:`j`, a minimum distance
    :math:`D_j`, can be found to this set of points.
    Renormalization is accomplished by replacing the
    level set value at this node :math:`\phi_j` 
    with :math:`D_j` multiplied by
    the sign of the previous value for the level set
    function. This method is fast and robust and
    reasonably accurate given sufficiently refined
    meshes using high order level set interpolation. 
    However, this method is prone to losing material if
    low order level set interpolation is employed.

Huygens_Constrained
    This method renormalizes the function in much the same way as the
    **Huygens** method, except it employs a Lagrange multiplier to enforce
    a global integrated constraint that requires the volume occupied by the
    “negative” phase to remain unchanged before and after renormalization. This
    requirement makes this method better at conserving mass. However, since it
    enforces a global constraint, it is possible that material might be moved
    nonphysically around the computational domain.

------------
Examples
------------

This is a sample renormalization method input card: 
::

	Level Set Renormalization Method = Huygens_Constrained

-------------------------
Technical Discussion
-------------------------

Renormalization is an operation particular to level set embedded interface tracking.
The level set function :math:`\phi` is usually specified in terms of a signed distance to the
interface. This type of function has very nice properties in terms of smoothness and a
unitary gradient magnitude in the vicinity of the interface. All of which are beneficial
in accurately integrating the function and applying interfacial physics such as surface
tension. The difficulty appears because of the velocity field :math:`\underline{u}` used to evolve the level
set function via the relation:

.. math::

   \frac{\partial \phi}{\partial t} + \underline{u} \cdot \nabla \phi = 0

There is nothing that requires that this velocity preserve the level set function as a
distance function during its evolution. As a result, large gradients in the level set
function might appear that would degrade the accuracy of both its time evolution and
the accuracy of the interfacial terms related to the level set function. To remedy this
problem, periodically the level set function must be reconstructed as a distance
function; this process is referred to as renormalization. The criteria for determining
when renormalization should occur is discussed under *Level Set Renormalization
Tolerance.*

