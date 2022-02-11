*********************************
**Number of Jacobian File Dumps**
*********************************

::

	Number of Jacobian File Dumps = <integer>

-----------------------
**Description / Usage**
-----------------------

This routine will dump a serial machine independent binary file out to disk containing
the Jacobian. The file is meant to be used by the auxiliary program, **checkGomaJac**, to compare two versions of the Jacobian. Ancillary data meant to enhance the printouts in
**checkGomaJac** are also output to the file. The card takes one mandatory integer
variable.

<integer>       
    If the integer is a positive number, n, then *Goma* will dump the first
    n Jacobians created (for any reason) to the current directory. If the
    integer is a negative value, -n, then *Goma* will dump a single Jacobian,
    the n’th Jacobian created, to the current directory.

The dumped files are named matrix.000, matrix.001, etc. Overwrites of files are
allowed to occur. The files themselves are written out using the XDR protocol layer
(easy, quick, and machine portable). The VBR format is used to write files out, even if
the internal format used by *Goma* is MSR. Thus, VBR and MSR formatted Jacobians
may be compared. Frontal Solver Jacobians are not compatible. The algorithm used is
also compatible with parallel jobs using *Goma*. In other words, the Jacobian file
dumped out for an 8 processor *Goma* run should be identical to the file dumped out by
a single processor run.

In order to use this feature, it is necessary to compile *Goma* with the MATRIX_DUMP
flag defined.

To compare two Jacobian files previously dumped out for compatibility, run
**checkGomaJac** offline:


::

	checkGomaJac   matrix1   matrix2

**checkGomaJac** will compare each entry in the row and column scaled matrices and
print out in an annotated format the entries containing the largest differences.

------------
**Examples**
------------

::

	Number of Jacobian File Dumps = 2

-------------------------
**Technical Discussion**
-------------------------

This capability has proven itself to be very useful in tracking changes to the Jacobian
due to differences in the machine architecture, number of processes, and due to changes
in the source code over time. The comparison is done using the standard RTOL, ATOL
logic found in ODE solvers. In other words, a weighting vector of the form,

**EQUATION** 

is created for each Jacobian entry, J\ :sub:`i`. Then, a determination of the difference between
J\ :sub:`i`and J\ :sub:`i` by the following formula:

**EQUATION**

w\ :sub:`i` is also used in the Jacobian column scalings, before the standard row sum scaling is
applied.

Internal Sandia users can find the auxiliary program, **checkGomaJac**, in the directory /
**home/goma/arch/linux/bin** on the Linux compute server, and in other ‘arch’
subdirectories for other platforms. External users should contact Goma support staff to
obtain the tool.




.. 
	TODO - I tried wriitng the equation out, but line 66 is definatly not how it is suppsoed to look. I could not for the life of me figure out how to do a superscript and cubscript on one letter in lines 68-69. Kris said there is a math thing to help me do this.
