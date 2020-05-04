# Goma 6.2
A Full-Newton Finite Element Program for Free and Moving Boundary Problems with Coupled Fluid/Solid Momentum, Energy, Mass, and Chemical Species Transport

For more information see the [Goma website](http://goma.github.io)

Changes in version 6.2

* Many bug fixes
* Stratimikos support
* Log conformation tensor stress model
* Hysing and Denner surface tension models for level set
* Suspension balance updates
* Updated SUPG for species
* Quadratic triangles
* And more...

Changes in version 6.1

* Automatic creation of brk files (see `docs/parallel_integration.md`)
* Bug fixes (mostly MPI related bugs)
* Epetra is now a supported matrix format (see `docs/epetra_matrix_format.md`)
* MUMPS is now supported through Amesos
* An experimental build script to build library dependencies is available under `scripts/`

## Build Instructions

### Get the Goma source code

Run the command

`git clone https://github.com/goma/goma`

### Get the prerequisites

Goma requires many third party libraries (TPLs, listed in [scripts/README.md](scripts/README.md)) which can be built automatically with the build script located in the scripts directory.

#### Build script requirements

The build-goma-dependencies.sh script relies on several packages readily available in many repositories.

For Ubuntu this will install the necessary packages to run the script:

`sudo apt-get install build-essential m4 zlib1g-dev libx11-dev gfortran`

For CentOS

`sudo yum install gcc-c++ gcc-gfortran m4 zlib-static libX11-devel`

The [scripts README](scripts/README.md) gives a more specific list of minimum requirements. Please ensure the script requirements are met before attempting to use it.

### Build Everything

To build third party libraries (TPLs):

`cd goma/scripts`

`./build-goma-dependencies.sh [-j NUM_PROCS] GOMA_TPL_INSTALL_DIRECTORY`

Where `-j NUM_PROCS` is optional to run the build with multiple processors and
`GOMA_TPL_INSTALL_DIRECTORY` is where you would like the third party libraries
to be installed.

Add the openmpi `bin` directory to your `PATH`, the trilinos `bin` directory to
your path, and the openmpi `lib` directory to your `LD_LIBRARY_PATH`
See [scripts/README.md](scripts/README.md) for more details on setting paths.

Finaly to build goma copy the `settings.mk-example` to `settings.mk` edit
`GOMA_LIBS` to point to your TPL directory, and make sure the trilinos and
openmpi paths/versions are correct.

Then to build goma run:

    make

The `goma` executable will be built in the `[path to goma]/bin` directory

Optionally `goma` can be installed with

    make install PREFIX="/path/to/install/"

default prefix if unspecified is `/usr/local`

## Environment Variables
OpenMPI should be added to the path and library path.
If the build script was used to build openmpi you would use something like:

    export LD_LIBRARY_PATH="/[path to gomalibs]/openmpi-4.0.2/lib:$LD_LIBRARY_PATH"
    export PATH="/[path to gomalibs]/openmpi-4.0.2/bin:$PATH"

SEACAS tools from Trilinos (e.g. `aprepro` and `blot`) should be added to
your path

    export PATH="/[path to gomalibs]/trilinos-12.18.1/bin:$PATH"  

For netcdf utilities such as `ncdump` add netcdf executables to your path

    export PATH="/[path to gomalibs]/netcdf-c-4.7.3/bin:$PATH"  

### Run the tutorial

To get started with Goma, use the following:

* [Tutorial instructions](http://goma.github.io/files/goma-beginners-tutorial.pdf)
* [Tutorial files tarball](http://goma.github.io/files/goma_beginners_tutorial.tar.gz)
