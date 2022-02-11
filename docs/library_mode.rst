======================================
Appendix 2: Using Goma in Library Mode
======================================

A new capability has been added to Goma which allows it to be linked with another finite element
program. This mode allows Goma to be compiled as a set of subroutines, which can be called
from another program. This will allow an external driver to use Goma and another code to solve a
problem which can be suitably decoupled - so that Goma solves some of the governing equations
and the other code solves the others. In this mode, there are pre-defined sets of variables or fields
which each code is responsible for assembling and passing on to the other code through a
common driver. This was designed so that Goma can be coupled with JAS3D; however, an
attempt has been made to develop a general capability which can be used with other codes as well.
This may minimize the need to implement new physics equations in Goma when other codes with
the desired routines are available.

A major component of this implementation within the Goma source is the addition of an alternate
version of the main program which is called “jas_main.c.” This version divides all tasks into three
subroutines:

* goma_init: Initializes code, parses input and broadcasts, do global array allocations
* goma_solve: Calls Goma transient solver
* goma_close: Cleans up

There is no "main" function, so there is no Goma executable as such. Instead, the source files are
compiled and assembled into the standard libraries libgoma.a and libgomau.a. These libraries,
along with those of the other program, are then used for linking the common driver, which may be
in a language other than C.

Communication between the codes is handled through four 1D arrays which are passed into and
out of Goma:

* xnv_in: Values of nodal variables imported into Goma.
* xev_in: Values of element variables imported into Goma.
* xsoln: Values of Goma solution variables exported from Goma.
* xpost: Values of Goma post processing variables exported from Goma.

These arrays can accommodate multiple variables, one right after the other: First x1[0..N], then
x2[0..N], and so forth. where N is the number of dofs of that variable in the problem, which is the
number of nodes except for xev_in where it is the number of elements. For Goma’s purposes,
these are considered external fields. For this reason, the MAX_EXTERNAL_FIELD in the Goma
makefile must be set high enough for the number of imported variables. The values passed in are loaded into the efv->ext_fld_ndl_val[] arrays, and used within Goma just as if the values were
read from an external file. Note that when mesh displacements are imported this way, Goma uses
a flag efv->ev_porous_decouple which must be set to TRUE - this signals Goma to anneal its
undisplaced mesh with the external displacements, so that the displaced nodal coordinates are
used for Jacobians, etc. without having to turn on Goma’s mesh equations.

For the time being, it is assumed that all external fields are used as nodal variables within Goma,
but may be element variables in the code that calculates them. Therefore, a routine has been added
to interpolate imported element variables to the nodes. This is a very naive linear interpolation,
but if this appears to be insufficient, a higher-order interpolation scheme can be easily
implemented later. The interpolated values still end up in the efv->ext_fld_ndl_val[] array(s).

The fields Goma will be importing must be specified in the input deck. To do this, use the same
"External Field" card as before, but instead of specifying a file name, place the string "IMPORT"
for nodal values or "IMPORT_EV" for element variables which must be interpolated. PLEASE
NOTE that the order of the External Field input cards will determine the order in which the values
must be loaded into the import arrays. Also, all nodal field cards must be placed before any
element field cards (otherwise, an EH will result).

There are two types of fields which can be exported from Goma: variables direct from the Goma
solution vector x[], and scalar post-processing variables. The convention for these is similar.
To specify a solution variable for export, add the following card below the last External Field
card:

**Export Field = <N>**

where <N> is the integer value assigned to the variable in the file rf_fem_const.h. There can be up
to MAX_EXTERNAL_FIELD of these variables to be exported, and they will be loaded into the
xsoln array in the order of the cards.

Exporting post-processing variables is a little more tricky. Initially, onlyscalar post-processing
fields have been enabled to be specified for export, in order to simplify the allocation process. The
export of vector fields (such as electric field) or tensor fields (such as stress) can be enabled in the
future as the need arises. To do this, go to the relevant card in the Post Processing Specifications
input section, and change the "yes" to "exp" to enable space for it to be allocated in the xpost
array. Note, however, that the order in which these fields will be stored in that array (when there
are two or more) is determined by the order in which they are processed in the function
load_nodal_tkn(), which may differ from the order of the cards in the input deck.

The way the post-processing export scheme works is as follows: There is a new array "x_pp" of
type double in solve_problem, which is allocated to size (NNODES *
MAX_EXTERNAL_FIELD) when the LIBRARY_MODE flag is defined or left NULL
otherwise. This array is passed into write_solution(), and in turn passed into
post_process_nodal(), so that once the post-processing fields are calculated, the requested values
can be saved there - otherwise they would simply be dumped into the Exodus file and erased from
memory. This is why there is now an extra argument to each of these functions. These saved
values are then loaded into xpost before x_pp itself is deallocated.

The four import/export arrays are intended to be allocated within the driver code, using
information obtained from parsing the Goma input deck and passed back to the driver upon
exiting goma_init(). This information consists of the number of fields to be stored in each array,
the number of elements, and the number of nodes. These are pointer arguments to the function
goma_init().

It is anticipated that Goma will be used in library mode to solve transient problems. Therefore, a
provision has been made for Goma to be called as a subroutine several times during a run, with a
start and end time passed in on each step. Goma may take one or several steps to reach the
requested end time on any given call. In any case, the actual end time is passed back to the driver,
with a warning if the requested time was not reached (e.g. due to Goma step failure), so that the
other code will know exactly how far to proceed to remain in sync with Goma. It is also possible
to have the other code precede Goma at each step. This is handled in the driver code, which passes
an argument to Goma indicating which code is called first.

To build the Goma libraries (libgoma.a and libgomau.a) in library mode, the makefile must be
modified as follows: Replace main.c and main.o with jas_main.c and jas_main.o in the
MAIN_SRC and MAIN_OBJ lists, and add the flag -DLIBRARY_MODE to the list of
DEFINES. This compiler flag activates many sections of code which were added in developing
this capability, and also invokes expanded argument lists for some functions which handle
communication data and arrays. Note that since there is no program “main” in jas_main.c, it will
not be possible to generate a stand-alone Goma executable in this way. This may cause an error
message on some platforms (even though the libraries are created successfully); to remedy this, it
is possible to create a new target in the makefile (e.g. “goma_jas”) in which the final command to
create the Goma executable is omitted.

Once the libraries for both Goma and the other program are built, then the driver code can be
compiled and linked with these libraries included to create a global executable. The driver
currently available for Goma is ANIMAS, which links to Goma and JAS3D. To obtain this driver
and specific build instructions, please contact Edward Wilkes (edwilke@sandia.gov).

