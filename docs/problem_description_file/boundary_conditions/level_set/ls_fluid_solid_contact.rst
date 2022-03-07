**************************
**LS_FLUID_SOLID_CONTACT**
**************************

::

	BC = LS_FLUID_SOLID_CONTACT LS <integer> <integer1>

-----------------------
**Description / Usage**
-----------------------

**(EMB/MOMENTUM)**

This boundary condition applies a fluid-solid stress balance at a level set interface that
is slaved to an overset mesh (see GT-026.3). It is applied as an “embedded” source
term on the fluid momentum equations at the zero level set contour. **NOTE**: *This
boundary condition has been deprecated in favor of the* **BAAIJENS_SOLID_FLUID**
*and* **BAAIJENS_FLUID_SOLID** *boundary conditions, as described in the memo*.

A description of the input parameters follows:

========================== ===============================================================
**LS_FLUID_SOLID_CONTACT** Name of the boundary condition.
**LS**                     This string is used to indicated that this is a “boundary”
                           condition is applied at an internal phase boundary defined
                           by the zero contour of the level set function.
<integer>                  An integer parameter than is permitted to take one of three
                           values -1, 0, or 1. Depending upon the choice of this
                           parameter the mass flux value is applied to the negative
                           phase, both phase, or the positive phase, respectively.
                           Details are given below.
<integer1>                 *Not used. Set to zero.* 
========================== ===============================================================

-------------------------
**Technical Discussion**
-------------------------

We discourage use of this experimental boundary condition.



--------------
**References**
--------------

No References. 




