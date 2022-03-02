***************
**SOLID_FLUID**
***************

::

	BC = SOLID_FLUID SS <bc_id> <integer1> <integer2> [float]

-----------------------
**Description / Usage**
-----------------------

**(PCC/VECTOR REALSOLID)**

The *SOLID_FLUID* condition performs the exact same task as the *FLUID_SOLID*
condition. The usage and example are also the same, so consult the discussion on that
card for further information.

At one time this condition applied the stress balance between solid and fluid phases in a
different fashion that proved not to be useful. To preserve backward compatibility, we
have kept this boundary condition around even though it invokes the exact same
function that the *FLUID_SOLID* boundary condition does.

Definitions of the input parameters are as follows:

=============== =============================================================
**SOLID_FLUID** Name of the boundary condition (<bc_name>).
**SS**          Type of boundary condition (<bc_type>), where **SS** denotes 
                side set in the EXODUS II database.
<bc_id>         The boundary flag identifier, an integer associated with
                <bc_type> that identifies the boundary location (side set in
                EXODUS II) in the problem domain.
<integer1>      Element block ID of solid phase from the EXODUS II
                database.
<integer2>      Element block ID of liquid phase from the EXODUS II
                database.
[float]         Scale factor for stress balance for non-dimensionalization.
                This parameter, which multiplies the liquid phase
                contribution of stress, is optional. The default is 1.0.
=============== =============================================================

------------
**Examples**
------------

See *FLUID_SOLID* description.

-------------------------
**Technical Discussion**
-------------------------

See *FLUID_SOLID* description.



--------------
**References**
--------------

No References.