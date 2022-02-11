==========================
Code Structure and I/O
==========================

Files for Data Input
########################

The *Goma* file I/O structure is diagrammed in Figure 2. Input to the program is divided into six
categories: (1) command-line options, (2) problem description file, (3) material files, (4) ASCII
continuation/restart file, (5) EXODUS II database file, and (6) sundry material property or
boundary condition table lookup files. *Goma* is basically set up to run in batch mode, i.e., no input
is required on the command line or after the run command is issued. There are, however, several
command-line switches which can be used to redirect I/O, control the level of I/O, and activate
debugging options.

.. figure:: /figures/002_goma_physics.png
	:align: center
	:width: 90%

	I/O structure for *Goma*. Dashed lines indicate that the files or commands are
	not required.

The *problem-description* file is by default called “input” but can be renamed with the -i switch on
the command line. A version of this file is also output as an “echo” file, viz. a prefix “echo”
prepended to the input file name. The echo file is used to verify input into goma, as it clearly
states all default settings for the input file and material files. . The input file itself contains the
general description of the problem and directions to *Goma* on how to solve it (see Chapter 4). The
file is split into thirteen sections: (1) File Specifications (Section 4.1) which directs I/O, (2)
General Specifications (Section 4.2), (3) Time Integration Specifications (Section 4.3), (4)
Continuation Specifications (Section 4.4), (5) Hunting Specifications (Section 4.5), (6)
Augmenting Condition Specification (Section 4.6), (7) Solver Specifications (Section 4.7), (8) Eigensolver Specifications (Section 4.8), (9) Geometry Specification (Section 4.9), (10)
Boundary Condition Specifications (Section 4.10), (11) Rotation Specifications (Section 4.11),
(12) Problem Description (Section 4.12), and (13) Post Processing Specifications (Section 4.13);
this latter section includes breakouts for fluxes and data (Section 4.14), particle traces (Section
4.15) and for volume-based integrals. The file format is described in detail in Chapter 4.
Incidentally, the structure of the data input routines is divided roughly along the same lines as the
input data file itself. 

The *material description* files (using the nomenclature “[material name].mat”) contain all
material property data and material property model and constitutive model specifications. The
names of these files are specified in the problem description file. The format of these files and the
available options are described in Chapter 5. Note that these files are also reproduced as output as
“echo” files, with all default settings specified.

The *ASCII continuation/restart files* (may have any name) contain an ASCII list of the solution
vector (values of field variables at nodes), which can be used as an initial guess for successive
runs of *Goma*. The names of these files are specified in the problem description file, but may be
changed with the -c (for input) or -s (for output) command-line options. These restart files are
“recyclable”, in the sense that output from one *Goma* simulation may be used as input to another
*Goma* simulation under certain restrictions.

The *EXODUS II database files* (may have any name but generally end in “.exoII”) contain a
description of the finite-element structure for the current problem. All EXODUS II files contain a
definition of the mesh, material blocks, and boundary sets. In the case of input EXODUS II files
created from mesh generator output, this is the sole content of the file. Output EXODUS II
database files contain a clone of the input EXODUS II mesh information and also contains the
nodal values of all field variables in the solution. The names of these files are specified in the
problem description file, but may be changed with the -ix (for input) or -ox (for output)
command-line options. The only EXODUS II file required when running *Goma* is the one
containing the current problem mesh. All others are either output for postprocessing or used to
supply auxiliary external fields (e.g. magnetic fields).

Command-Line Arguments
##########################

*Goma* can be run using only the input files (all four listed above) to describe the problem and to
direct the input and output; in this case *Goma* is run using the command “goma” without any
arguments. However, command-line arguments offer additional flexibility for redirecting input or
output and for adjusting common run-time parameters. The general command line for running
*Goma* is:

:: 

    $ goma [-nd] [-se fn] [-so fn] [-i fn] [-c fn] [-s fn] [-ix fn] [-ox fn] [-d int] 
           [-n int] [-r dbl] [-a args] [-restart fn] [-h] [-ts dbl] [-te dbl] [-cb dbl] 			
           [-ce dbl] [-cd dbl] [-cn int] [-cmin dbl] [-cmax dbl] [-cm int] [-ct int] [-c_bc int]
           [-c_df int] [-c_mn int] [-c_mp int] [-bc_list] [-v]

Here *fn* denotes “file name”, *int* denotes “integer”, *dbl* denotes “float or double” and *args* denotes
multiple sub-options or file names. The input line is parsed into options, which are preceded by a
single hyphen (-) and arguments, which normally are *fn*, *int*, or *dbl* not preceded by a hyphen. The
default, if no options are specified, is the input option (e.g. “goma input.alt” is the same as
“goma -i input.alt”). The following is a list of the command-line options and their
descriptions (two ways are shown to specify each option, an abbreviated and a verbose form).

-a args, -aprepro args    
                          Preprocess input files through the APREPRO preprocessor 
                          [with args as arguments to APREPRO] before reading into 
                          *Goma*. With this option, *Goma* performs a UNIX
                          system() call to run APREPRO which will preprocess
                          the input file and the material data files. The
                          APREPRO input file is preprocessed from “input” or
                          the filenamespecified by the -input option and
                          written to “tmp.input”. Likewise, the material
                          data files are preprocessed from “[material
                          name].mat” to “tmp.[material name].mat”. After the
                          “-a” on the command line, options for APREPRO are
                          preceded by two hyphens (--). For example, the
                          command line “goma -i *input.pre* -a CONSTANT1=0.2
                          --vd” will preprocess “input.pre” and the material
                          data files specified in *input.pre* using APREPRO,
                          and will pass the argument -*vd* (which prints
                          version number and values of all variables to the
                          screen) and CONSTANT1=0.2 (which sets the variable
                          CONSTANT1 equal to 0.2 for preprocessing) to
                          APREPRO; the preprocessed files will be
                          “*tmp.input*” and “*tmp*.[material name].*mat*”.)

-c fn -contin fn          
                          Change the name of the ASCII continuation/restart input 
                          file (specified in Problem-Description File) to *fn*,
                          (e.g. “goma -c *old.soln.dat*” uses the file
                          “*old.soln.dat*” as the ASCII input file). Note that
                          this option has no effect if the initial guess is not
                          read from the ASCII file, i.e. unless “*Initial Guess
                          = read*” is specified in the input file.

-d int -debug int         
                          Change the debug flag to *int*. This option is convenient 
                          when debugging and the user wants to see more output
                          from *Goma*. (e.g. “goma - d -2” will run *Goma* with
                          the Debug_Flag set to -2). Higher values generally
                          produce more output.

-h -help                  
                          Prints a helpful message with brief descriptions of 
                          these command line options.

-i fn -input fn           
                          Redirect *Goma* to read the problem description file from 
                          *fn*. The normal default option is to read from
                          a file named “*input*”.

-ix fn -inexoII fn        
                          Redirect *Goma* to read the input EXODUS II database
                          file (often called “*in.exoII*”) from *fn*.

-brk fn
                          Specify a Brk file, this overrides setting Brk File in
                          Problem Description File

-n int                    
                          Change the maximum number of Newton iterations to *int*. This is
                          especially convenient for setting the number of
                          iterations to zero so that *Goma* just runs the
                          post-processor on the set of input data.

-nd -nodisplay            
                          Do not display the run-time information on the
                          screen. With this option, *Goma* sends the stdout and
                          stderr output to temporary files that are removed at
                          the end of the run. This command takes no arguments.

-ox fn -outexoII fn       
                          Redirect *Goma* to write the output EXODUS II file
                          (often called “*out.exoII*”) to *fn*.

-r dbl relax dbl          
                          Change the value of the Newton relaxation parameter
                          to *dbl*. This is convenient if a few Newton steps
                          with relaxation are desired before using full Newton.
                          (e.g. “goma -r 0.1” will use Newton’s method with
                          updates one-tenth of the normal value.

-s fn -soln fn            
                          Redirect Goma to write the output ASCII file
                          (normally called “*soln.dat*”) to *fn*.

-se fn -stderr fn         
                          Redirect the standard error from *Goma* to *fn*. This
                          output is comprised of more urgent diagnostic error
                          and timing messages.

-so fn -stdout fn         
                          Redirect the standard output from *Goma* to *fn*.
                          This output is comprised of less urgent informational
                          messages.

-ts dbl                   
                          Start time of simulation.

-te dbl                   
                          End time of simulation

-cb dbl                   
                          Continuation: Start value (see Gates et al., SAND2000-2465)

-ce dbl                   
                          Continuation: Final value (see Gates et al., SAND2000-2465)

-cd dbl                   
                          Continuation: Path step, ds (see Gates et al., SAND2000-2465)

-cn dbl                   
                          Continuation: Max number of path steps (see Gates et al., 2000)

-cm int                   
                          Continuation: Method (see Gates et al., 2000)

-ct int                   
                          Continuation: Type (see Gates et al., 2000)

-c_bc int                 
                          Continuation: Boundary condition ID (see Gates et al., 2000)

-c_df int                 
                          Continuation: BC Data Float ID (see Gates et al., 2000)

-c_mn int                 
                          Continuation: Material ID (see Gates et al., 2000)

-c_mp int                 
                          Continuation: Method property ID (see Gates et al, 2000)

-bc_list                  
                          Continuation: Method property ID (see Gates et al, 2000)

-v -version               
                          Output goma version

.. NOTE:: To get the most up-to-date list, simple issue the* “goma -h” *command
   at the command line. Also note that the continuation input parameters are
   explained in the Advanced Capabilities Manual (Gates et al. 2000 or newer
   version).*

The primary purpose of the command-line options is to allow the user an easy way to redirect the
input and output of Goma or to quickly change problem specifications. Most of the options are
overrides of information in the problem description file, so in some cases it may be easier to edit the problem description file than to use command-line arguments.

.. 
	TODO - In this document when it starts with "**-a** *args* **-aprepro** *args*" needs help being formatted. The table is smashing the word too much, so the idea is unclear. It should be 3 column to help distinguish, but the width of the column need to be figured out. 
