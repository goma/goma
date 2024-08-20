# Goma 
A Full-Newton Finite Element Program for Free and Moving Boundary Problems with Coupled Fluid/Solid Momentum, Energy, Mass, and Chemical Species Transport

For more information see the [Goma website](https://www.gomafem.com)

## Documentation

Most of the documentation can be found at [https://www.gomafem.com/documentation.html](https://www.gomafem.com/documentation.html)

## License

See [LICENSE](LICENSE) file. 
and are noted at the top of the cmake file.

### Third party library licenses

#### CMake modules 

Some cmake modules under `cmake/` were modified from the Eigen library
and are noted at the top of the cmake file.

See licenses at https://gitlab.com/libeigen/eigen

FindMETIS.cmake

* @copyright (c) 2009-2014 The University of Tennessee and The University
  of Tennessee Research Foundation. All rights reserved.
* @copyright (c) 2012-2014 Inria. All rights reserved.
* @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.


FindUMFPACK.cmake


#### nanoflann is included under the BSD license, please see `nanoflann.hpp`

 * Copyright 2008-2009  Marius Muja (mariusm@cs.ubc.ca). All rights reserved.
 * Copyright 2008-2009  David G. Lowe (lowe@cs.ubc.ca). All rights reserved.
 * Copyright 2011-2022  Jose Luis Blanco (joseluisblancoc@gmail.com).

## Major Changes

See [CHANGES.md](CHANGES.md)

## Build Instructions

See [BUILD.md](BUILD.md)

## Spack package

The Spack package manager [https://spack.io](https://spack.io/) can be used to install 
Goma and all of Goma's third party libraries

Currently available on the `develop` branch of spack.

Example for a bash-like shell:

    git clone https://github.com/spack/spack.git
    . spack/share/spack/setup-env.sh
    spack install goma

For more information on build options see:

    spack info goma

For more information on using spack see the [spack documentation](https://spack.readthedocs.io/en/latest/).


## Third party libraries

- Metis 5.1.0 (Optional)
- SEACAS 2022-01-27 (Required: Exodus and Aprepro)
- BLAS/LAPACK (Configured through Trilinos)
- Trilinos matrix solvers 13.0.1 and up (Required: AztecOO, Amesos, Epetra, TPL LAPACK; Optional: Stratimikos [with Teko, Ifpack, Belos, Tpetra])
- PETSc matrix solvers (KSP, PC)
- MUMPS 5.4.0 (through Trilinos or PETSc only)
- Superlu_dist 7.2.0 (through Trilinos or PETSc only, Trilinos requires parmetis build)
- UMFPACK, SuiteSparse 5.10.1 (Optional)
- ARPACK/arpack-ng 3.8.0 (Optional)
- sparse 1.4b (Optional)
- Catch2 (Optional testing)

### Run the tutorial

To get started with Goma, use the following:

* [Tutorial instructions](https://docs.gomafem.com/files/goma-beginners-tutorial.pdf)
* [Tutorial files tarball](https://docs.gomafem.com/files/goma_beginners_tutorial.tar.gz)
