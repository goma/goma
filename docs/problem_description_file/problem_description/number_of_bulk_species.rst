**************************
**Number of Bulk Species**
**************************

::

	Number of bulk species = <integer>

-----------------------
**Description / Usage**
-----------------------

This card is required for each material section in the *Problem Description File*. It is
used to specify the number of species in a phase. The word, bulk, here, refers to its
being distributed throughout the domain, not just at a surface. All loops over property
evaluations use this value to specify the length of the loop. The single input parameter
is defined as:

========= ====================================================================
<integer> The number of species. If the value of <integer> is 0, then
          no species equations are solved for.
========= ====================================================================

In the absence of any further cards specifying the number of species equations, the
number of species equations is set equal to the integer value supplied by this card, and
there is an implied additional species, i.e., the solute, which is not part of species
loops, but which fills out the specification of the phase.

------------
**Examples**
------------

Following is a sample card:
::

   Number of bulk species = 1

-------------------------
**Technical Discussion**
-------------------------

Unfortunately, in the past, this card has specified the number of species equations
instead of the number of species, as its name would imply! Now, the preferred
treatment is to specify unequivocally both the number of bulk species and the number
of bulk species equations using two separate input cards. If the two values are the same,
then the system is semantically referred to as being “dilute” (even though it might not
be!), and there is an inferred solute which is not part of the loop over species unknowns
in property evaluations or even in the specification of properties in the .mat file. If the number of species is one greater than the number of species equations, then the system
is deemed “nondilute” and the length of loops over property evaluations is one greater
than the number of species equations. For nondilute systems, an equation of state must
be implicitly used within *Goma* to solve for the value of the species unknown variable
for the last species.



--------------
**References**
--------------

No References.