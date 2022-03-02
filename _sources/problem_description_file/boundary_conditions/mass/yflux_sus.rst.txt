*************
**YFLUX_SUS**
*************

::

	BC = YFLUX_SUS SS <bc_id> <integer>

-----------------------
**Description / Usage**
-----------------------

**(WIC/MASS)**

This boundary defines a flux of suspension particles at an interface. Definitions of the input parameters are as follows:

============== =================================================================
**YFLUX_SUS**  Name of the boundary condition (<bc_name>).
**SS**         Type of boundary condition (<bc_type>), where **SS** denotes
               side set in the EXODUS II database.
<bc_id>        The boundary flag identifier, an integer associated with
               <bc_type> that identifies the boundary location (side set in
               EXODUS II) in the problem domain.
<integer>      Species id; the species number for suspension particles.
============== =================================================================

------------
**Examples**
------------

The following is a sample input card:
::

   BC = YFLUX_SUS SS 1   0

-------------------------
**Technical Discussion**
-------------------------

This condition is only used in conjunction with the *SUSPENSION* liquid constitutive
models, *HYDRODYNAMIC* diffusivity model, and *SUSPENSION* or
*SUSPENSION_PM* density models. A theoretical outflux condition associated with
suspension particles leaving the domain is tied to the Phillips diffusive-flux model.
Please refer to discussions on *HYDRODYNAMIC* diffusivity to gain more
understanding of the suspension flux model.



