************
**FRICTION**
************

::

	BC = FRICTION SS <bc_id> <float1> [integer1]>

-----------------------
**Description / Usage**
-----------------------

**(WIC/VECTOR MESH)**

This boundary condition card applies a force per unit area (traction) on a Lagrangian
mesh region. The force per unit area is applied according to Coulomb’s friction law
over the boundary delineated by the side set ID. The applied traction is of course a
vector. The vector traction is defined in normal-tangent vector basis. Definitions of the
input parameters are as follows:

============ ========================================================
**FRICTION** Name of the boundary condition (<bc_name>).
**SS**       Type of boundary condition (<bc_type>), where **SS**
             denotes side set.
<bc_id>      The boundary flag identifier, or a side set number which
             is an integer that identifies the boundary location (side
             set in EXODUS II) in the problem domain.
<float1>     :math:`\mu`, Coulombic coefficient of friction.
[integer1]   optional specification of the element block id to which
             this condition will be applied.
============ ========================================================

This card actually applies a traction that is then naturally integrated over the entire side
set of elements.

------------
**Examples**
------------

Following is a sample card:
::

     BC = FRICTION SS 10 0.1 2

-------------------------
**Technical Discussion**
-------------------------

**Important note**: this boundary condition can only be applied to *LAGRANGIAN,
DYNAMIC_LAGRANGIAN* or *ARBITRARY* mesh motion types (cf. *Mesh Motion* card).
For real-solid mesh motion types, refer to *FRICTION_RS*.

This condition should be utilized in conjunction with a rotated condition such as
SPLINE in order to apply a tangential force which is proportional to the normal force;

.. figure:: /figures/069_goma_physics.png
	:align: center
	:width: 90%

where :math:`\mu` is the coefficient of friction and :math:`\underline{v}` is the velocity of the convected solid. Note
that the direction of the frictional force is determined by the velocity direction.




.. 
	TODO - The image in line 56 needs to be replaced with the correct equation.