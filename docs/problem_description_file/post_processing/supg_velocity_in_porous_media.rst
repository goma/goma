*********************************
**SUPG Velocity in Porous Media**
*********************************

::

   SUPG Velocity in porous media = {yes | no}

-----------------------
**Description / Usage**
-----------------------

Used to specify use of effective velocities in SUPG formulations for porous media. It is
written to the output EXODUS II file as nodal variable **U_supg_porous**.

The permissible values for this postprocessing option are:

============= ================================================================
**yes**       Calculate and write the effective velocity components as a
              postprocessing variable to the output EXODUSII file.
**no**        Do not calculate the effective velocity components.
============= ================================================================

------------
**Examples**
------------

This is a sample input card to activate calculation of SUPG Velocity:
::

   SUPG Velocity in porous media = yes

-------------------------
**Technical Discussion**
-------------------------

No Discussion.



--------------
**References**
--------------

GTM-029.0: SUPG Formulation for the Porous Flow Equations in Goma, H. K.
Moffat, August 2001 (DRAFT).