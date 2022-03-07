**************************
**SURFACE_ELECTRIC_FIELD**
**************************

::

	BC = SURFACE_ELECTRIC_FIELD SS <bc_id> <integer> <integer> <integer>

-----------------------
**Description / Usage**
-----------------------

**(WSG/SURFACE CHARGE)**

This boundary condition card is used to apply a part of the shell surface charge
equation which includes the electric field, the negative gradient of the potential
variable which is applied on a neighboring bulk block. It is actually an integral part of
the surface charge equation. Definitions of the input parameters are as follows:

========================== ====================================================
**SURFACE_ELECTRIC_FIELD** Name of the boundary condition (<bc_name>).
**SS**                     Type of boundary condition (<bc_type>), where **SS**
                           denotes side set in the EXODUS II database.
<bc_id>                    The boundary flag identifier, an integer associated 
                           with <bc_type> that identifies the boundary location (side set in EXODUS II) in the problem domain. This boundary must coincide with the shell element block on which the surface charge equation is applied.
<integer>                  Bulk element block ID (from ExodusII database) for
                           neighboring bulk element block on which the potential equation is applied.
<integer>                  Shell element block ID (from ExodusII database) for
                           shell block on which the surface charge equation is
                           applied.
========================== ====================================================

This boundary condition is currently inoperative...

------------
**Examples**
------------

For a system consisting of a solid material (element block ID 1) with a conducting shell surface (element block ID 2) whose location coincides with side set 20, the following is a sample usage:
::

   BC = SURFACE_ELECTRIC_FIELD SS 20 1 2

-------------------------
**Technical Discussion**
-------------------------

This is a special type of boundary condition, WEAK_SHELL_GRAD, which is a
portion of a shell equation which involves spatial gradients of bulk variables. Since the
values of bulk variable gradients depend on all of the degrees of freedom of that
variable in the bulk element, and sensitivities to the off-shell degrees of freedom must
be applied, a portion of the equation must be evaluated from the bulk side. This is done
in Goma by means of a WEAK_SHELL_GRAD boundary condition which evaluates
these terms and all bulk sensitivities from the bulk side, the saves these values for later
recall when the rest of the surface charge equation is assembled.

In this case, the term n • :math:`\Delta` *V*
and its potential sensitivities are evaluated within the bulk
element for inclusion in the surface charge balance along the shell surface.. I



--------------
**References**
--------------

No References.