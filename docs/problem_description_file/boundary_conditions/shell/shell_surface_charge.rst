************************
**SHELL_SURFACE_CHARGE**
************************

::

	BC = SHELL_SURFACE_CHARGE SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(WIC/POTENTIAL)**

This boundary condition card is used to add to the potential equation the surface charge
term at a shell surface. Definitions of the input parameters are as follows:

======================== ========================================================
**SHELL_SURFACE_CHARGE** Name of the boundary condition (<bc_name>).
**SS**                   Type of boundary condition (<bc_type>), where **SS**
                         denotes side set in the EXODUS II database.
<bc_id>                  The boundary flag identifier, an integer associated with
                         <bc_type> that identifies the boundary location (side set
                         in EXODUS II) in the problem domain. This boundary
                         must coincide with the shell element block on which the
                         surface charge equation is applied.
<integer>                Integer value indicating the bulk element block ID from
                         which to apply the boundary condition (not currently
                         implemented).
======================== ========================================================

This boundary condition is currently inoperative.

------------
**Examples**
------------

For a system consisting of a solid material (element block ID 1) with a conducting shell
surface (element block ID 2) whose location coincides with side set 20, the following is
a sample usage:
::

   BC = SHELL_SURFACE_CHARGE SS 20 1

-------------------------
**Technical Discussion**
-------------------------

This boundary condition was originally developed to allow for fluid slip near a
dynamic contact line, a necessary condition for dynamic wetting line motion when the
contact angle is not 180 degrees (viz. rolling motion condition). The slippage
mechanism was deployed through the use of Navier’s slip condition, which basically
goes as

.. figure:: /figures/245_goma_physics.png
	:align: center
	:width: 90%

where **E** is the electric field vector, the superscripts (*o*) and (*i*) denote the outer and inner phases, *n* is a unit normal pointing into the outer phase,
:math:`\varepsilon` is the electrical permittivity, 
E = –:math:`\Delta` V is the electric field and *V* is the voltage or electric potential.




.. TODO -Line 55 has a photo hat needs to be replaces with an equation.
