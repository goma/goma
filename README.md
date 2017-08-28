# Goma 6.1
![](http://129.24.4.70:8080/buildStatus/icon?job=GomaMaster "Build Status")

A Full-Newton Finite Element Program for Free and Moving Boundary Problems with Coupled Fluid/Solid Momentum, Energy, Mass, and Chemical Species Transport

For more information see the [Goma website](http://goma.github.io)

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

The build-goma-dep-trilinos-12.sh script relies on several packages readily available in many repositories.

For Ubuntu this will install the necessary packages to run the script:

`sudo apt-get install build-essential m4 zlib1g-dev libx11-dev gfortran`

For CentOS

`sudo yum install gcc-c++ gcc-gfortran m4 zlib-static libX11-devel`

The [scripts README](scripts/README.md) gives a more specific list of minimum requirements. Please ensure the script requirements are met before attempting to use it.

### Build Everything

Simply use these commands:

`cd goma/scripts`

`./easy-goma-builder.sh`

and then follow the on-screen instructions to build Goma.

If this does not work, or you need something more specific (like a debugging version of Goma, or one that works without c++11 support), please follow the [manual build instructions](manualbuild.md)

### Run the tutorial

To get started with Goma, use the following:

* [Tutorial instructions](http://goma.github.io/files/goma-beginners-tutorial.pdf)
* [Tutorial files tarball](http://goma.github.io/files/goma_beginners_tutorial.tar.gz)
