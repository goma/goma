**********************
**DISCONTINUOUS_VELO**
**********************

::

	BC = DISCONTINUOUS_VELO SS <bc_id> <char_string> <integer1> <integer2>

-----------------------
**Description / Usage**
-----------------------

**(SIC/MOMENTUM)**

This boundary condition card, used to set the normal component of mass averaged
velocity at an interface, specifies that the net flux of the last component in a nondilute
mixture across an internal interface, is equal to zero. The condition only applies to
interphase mass, heat, and momentum transfer problems applied to nondilute material
phases with discontinuous (or multivalued) variables at an interface, and it must be
invoked on fields that employ the **Q1_D** or **Q2_D** interpolation functions to “tie”
together or constrain the extra degrees of freedom at the interface in question.

Definitions of the input parameters are as follows:

+---------------------+--------------------------------------------------------------+
|**DISCONTINOUS_VELO**| Name of the boundary condition (<bc_name>).                  |
+---------------------+--------------------------------------------------------------+
|**SS**               | Type of boundary condition (<bc_type>), where **SS**         |
|                     | denotes side set in the EXODUS II database.                  |
+---------------------+--------------------------------------------------------------+
|<bc_id>              | The boundary flag identifier, an integer associated with     |
|                     | <bc_type> that identifies the boundary location (side set    |
|                     | in EXODUS II) in the problem domain.                         |
+---------------------+--------------------------------------------------------------+
|<char_string>        | A character string identifiying the condition to be          |
|                     | applied on the liquid phase relative to the gas phase.       |
|                     |                                                              |
|                     |   * **EVAPORATION**                                          |
|                     |   * **DISSOLUTION** - not currently valid.                   |
|                     |                                                              |
|                     | Note, this parameter replaces the boundary condition         |
|                     | *EVAPORATION_VELO*.                                          |
+---------------------+--------------------------------------------------------------+
|<integer1>           | Element block id of liquid or high density phase.            |
+---------------------+--------------------------------------------------------------+
|<integer2>           | Element block id of gas or low density phase.                |
+---------------------+--------------------------------------------------------------+

------------
**Examples**
------------

Following is a sample input card that applies this BC on the block 1 side of side set 7,
the liquid side; the block 2 side is the gas side.
::

     BC = DISCONTINUOUS_VELO SS 7 EVAPORATION 1 2

-------------------------
**Technical Discussion**
-------------------------

The *DISCONTINUOUS_VELO* boundary condition applies the following equation:

.. figure:: /figures/094_goma_physics.png
	:align: center
	:width: 90%

It specifies the diffusive flux of the last species in the mechanism, i.e., the one for
which no explicit continuity equation exists, to be equal to zero. This is done via a
strong integral condition applied to one side of the interface, the “+” side of the
interface. This boundary condition, combined with the *KINEMATIC_SPECIES* and
*KINEMATIC_DISC* boundary conditions, implies that the diffusive flux of the last
species on both sides of the boundary is equal to zero.

The *DISCONTINUOUS_VELO* boundary condition requires an evaluation of the
derivative of the species mass fraction at the interface. Thus, the mesh convergence
properties of the algorithm are reduced to O( *h* ). Also, discretization error must interfere
with the total mass balance across a phase, since the expression for 
:math:`j_i^+` is substituted for
in some places, the *YFLUX_SPECIES* boundary condition, but used in the
*DISCONTINUOUS_VELO* boundary condition.




.. TODO - Line 65 contains a photo that needs to be exchanged for the equation.