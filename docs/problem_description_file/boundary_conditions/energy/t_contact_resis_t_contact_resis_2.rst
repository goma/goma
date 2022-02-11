**************************************
**T_CONTACT_RESIS, T_CONTACT_RESIS_2**
**************************************

::

	BC = T_CONTACT_RESIS SS <bc_id> <int1> <int2> <float1>

::

	BC = T_CONTACT_RESIS_2 SS <bc_id> <int2> <int1> <float1>

-----------------------
**Description / Usage**
-----------------------

**(WIC/ENERGY)**

This boundary condition set is used to specify a thermal contact resistance at an
interface between two mesh regions defined by a side set. Please see special usage
notes below regarding proper side-set specification and the reasons that both BC cards
are required for an internal interface. NOTE that the temperature field MUST be
interpolated with the discontinous versions of Q1 or Q2, viz. Q1_D and Q2_D.
Definitions of the input parameters are as follows:

=================== ==================================================================
**T_CONTACT_RESIS** Name of the boundary condition (<bc_name>).
**SS**              Type of boundary condition (<bc_type>), where **SS** denotes
                    side set in the EXODUS II database. Note this side set
                    MUST contain elements on both sides of the interface.
<bc_id>             The boundary flag identifier, an integer associated with
                    <bc_type> that identifies the boundary location (side set in
                    EXODUS II) in the problem domain.
<int1>              Material/region ID associated with first material.
<int2>              Material/region ID associated with second material. Note
                    that these material IDs are reversed on the second BC.
<float1>            Value contact resistance in units of thermal conductivity
                    divided by length.
=================== ==================================================================

------------
**Examples**
------------

The following is a sample card:
::

   BC = T_CONTACT_RESIS SS 3 1 2 10.0

::

   BC = T_CONTACT_RESIS_2 SS 3 2 1 10.0

Note that both boundary condition cards are required at an internal interface. In this
case the interface divides mesh/material ID 1 and 2. Note also how these material IDs
are reversed on the second card. These conditions apply a thermal contact resistance of
10. (units of thermal conductivity divided by length) at the interface defined by SS 3.

-------------------------
**Technical Discussion**
-------------------------

The mathematical form of the boundary condition.

.. figure:: /figures/134_goma_physics.png
	:align: center
	:width: 90%

The flux into the interface from material “a” is equivalent to that into material “b”, both equal to the temperature jump across the interface times the contact resistance 
:math:`R^{-1}`.

The side set to which this boundary condition is applied must contain elements on both
sides of the interface. Look up any special commands in your mesh generator to make
sure this occurs. In CUBIT, for example, you have to add “wrt volume 1 2” like
qualifiers on the side set command. The reason for the “double application” of this
condition is to pick up the all the terms from both sides of the interface with the proper sign. The nodes at the interface have two temperatures, one from each side, and so two weak form applications of this equation are required, one from each side.




.. TODO -Line 65 has a picture that needs to be changed out with the correct equation.