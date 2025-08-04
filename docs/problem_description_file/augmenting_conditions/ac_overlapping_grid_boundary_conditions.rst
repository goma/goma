******************************************
AC (Overlapping Grid Boundary Conditions)
******************************************

::

    AC = OV <SSID> <integer1> <integer2><integer3>

-----------------------
Description / Usage
-----------------------

This type of augmenting condition is used to invoke the overlapping grid algorithm for 
fluid-structure interaction problems, which applies a kinematic constraint along either 
the moving part of the solid boundary or the zero contour of a phase function field 
within the fluid which follows this boundary, and enables sensitivities from both phases 
to be included. When this card is provided in conjunction with an appropriate set of 
interfacial boundary conditions (see below), this AC is replaced with m new 
augmenting constraints -- two per element side along the side set specified by <SSID>.

This algorithm works with either of the following two sets of interfacial boundary 
conditions:

- LS_NO_SLIP, SOLID_FLUID_CONTACT, FLUID_SOLID_CONTACT

or

- LAGRANGE_NO_SLIP, OVERSET(BAAIJENS)_SOLID_FLUID, 
  OVERSET(BAAIJENS)_FLUID_SOLID

The primary difference between these two sets is that the former imposes the kinematic 
condition from the fluid side and the latter from the solid side. The same physical 
equations are applied in either case.

A description of the card syntax follows:

OV
    A mandatory string indicating that the augmenting 
    condition is being used to invoke the overlapping grid 
    algorithm.

<SSID>
    The SS index of the moving solid boundary, taken directly
    from the ExodusII input file.

<integer1>
    The element block index for the solid phase, taken directly 
    from the ExodusII input file.

<integer2>
    The element block index for the fluid phase, taken directly 
    from the ExodusII input file.

<integer3>
    Location of the Lagrange multiplier unknowns. For now, 
    this entry must be zero, which indicates that the unknowns 
    are stored in the augmenting conditions.

------------
Examples
------------

For a case where a fluid is in element block 1 and the solid is in element block 2 with 
sideset 10 defined along the moving boundary, the appropriate augmenting condition 
card of type OV is:

::

    AC = OV 10 2 1 0

Note that the first three integers come directly from the input mesh file and thus do not 
depend on which of the above boundary condition sets is chosen. It is only necessary 
that those BC's which are applied from the solid side have the same side set ID 
specified.

-------------------------
Technical Discussion
-------------------------

- If it is necessary to specify any other type of augmenting condition in the same 
  problem, this card must be the last one in the list.

- When this card is present, a check will be done for the presence of any of the above 
  interfacial boundary conditions. An error will occur if none are specified.

- There is a routine which counts the number of element sides in side set <SSID> 
  and allocate the necessary number of augmenting constraints to be used in the 
  program (This changes nAC). Goma must be compiled with a sufficient value for 
  the constant MAX_NGV, which must be at least nAC+5. The default value of this 
  constant is set to 10 in the file rf_io_const.h; this value can either be changed there 
  or overridden by specifying "-DMAX_NGV=<#>" in the DEFINES section of the 
  makefile.

--------------
References
--------------

Baaijens, F. P. T. "A fictitious domain/mortar element method for fluid-structure 
interaction," Int. J. Numer. Meth. Fluids. 35, 2001, 743-761.