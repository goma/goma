# Goma 
A Full-Newton Finite Element Program for Free and Moving Boundary Problems with Coupled Fluid/Solid Momentum, Energy, Mass, and Chemical Species Transport

For more information see the [Goma website](https://docs.gomafem.com)

## Documentation

Most of the documentation can be found at [docs.gomafem.com](https://docs.gomafem.com)

## License

See [LICENSE](LICENSE) file. Some cmake modules under `cmake/` were modified from the Eigen library
and are noted at the top of the cmake file.

## Major Changes

[Major changes for current version](CHANGES.md)

## Build Instructions

[Build instructions](BUILD.md)

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
