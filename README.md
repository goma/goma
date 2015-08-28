# Goma 6.1

A Full-Newton Finite Element Program for Free and Moving Boundary Problems with Coupled Fluid/Solid Momentum, Energy, Mass, and Chemical Species Transport

For more information see the [Goma website](http://goma.github.io)

Changes in version 6.1

* Automatic creation of brk files (see `docs/parallel_integration.md`)
* Bug fixes (mostly MPI related bugs)
* Epetra is now a supported matrix format (see `docs/epetra_matrix_format.md`)
* MUMPS is now supported through Amesos
* An experimental build script to build library dependencies is available under `scripts/`
* Goma now uses a settings.mk file to manage Makefile settings
