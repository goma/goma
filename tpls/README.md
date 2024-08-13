# Third party libraries

These are third party libraries with separate license files.

Please look at comments or license files for the separate license

# Script to install all tpls

Third party libraries required for Goma can be installed using the `install-tpls.py` script
available in this directory.

Examples:

Build with a nonstandard compiler:

    ./install-tpls.py --cc=gcc-13 --cxx=g++-13 --fc=gfortran-13 /path/to/install

Build with a prebuilt openmpi with 4 CPUs:

    ./install-tpls.py -j 4 --cc=mpicc --cxx=mpicxx --fc==mpifort --openmpi-dir=$MPI_HOME /path/to/install

Build with static libraries, shared is set by default:

    ./install-tpls.py --build-static /path/to/install

Reuse downloads:

    ./install-tpls.py --download-dir=/path/to/downloads /path/to/install

Default locations:

Downloads are downloaded to `/path/to/install/downloads`, and sources/builds are found in `/path/to/install/sources`
feel free to remove after installation.

## Usage

Once installed a `config.sh` and `config.fish` file are generated in `/path/to/install` and can be sourced
in your shell of choice and used to compile and run Goma.


Usage:
    
    usage: install-tpls.py [-h] [--cc CC] [--cxx CXX] [--fc FC] [--download-dir DOWNLOAD_DIR] [--extract-dir EXTRACT_DIR] [--build-shared] [--build-static] [-j JOBS] [--enable-parmetis]
                           [--disable-parmetis] [--cmake-dir CMAKE_DIR] [--openmpi-dir OPENMPI_DIR] [--hdf5-dir HDF5_DIR] [--pnetcdf-dir PNETCDF_DIR] [--netcdf-dir NETCDF_DIR] [--fmt-dir FMT_DIR]
                           [--seacas-dir SEACAS_DIR] [--bison-dir BISON_DIR] [--flex-dir FLEX_DIR] [--openblas-dir OPENBLAS_DIR] [--metis-dir METIS_DIR] [--parmetis-dir PARMETIS_DIR]
                           [--scotch-dir SCOTCH_DIR] [--arpack_ng-dir ARPACK_NG_DIR] [--scalapack-dir SCALAPACK_DIR] [--mumps-dir MUMPS_DIR] [--superlu_dist-dir SUPERLU_DIST_DIR]
                           [--suitesparse-dir SUITESPARSE_DIR] [--trilinos-dir TRILINOS_DIR] [--petsc-dir PETSC_DIR] [--petsc_complex-dir PETSC_COMPLEX_DIR] [--omega_h-dir OMEGA_H_DIR]
                           [--sparse-dir SPARSE_DIR]
                           INSTALL_DIR
    
    Third party library installer for the finite element code Goma
    
    positional arguments:
      INSTALL_DIR           Install location of TPLs
    
    options:
      -h, --help            show this help message and exit
      --cc CC               C compiler to use, default is CC environment variable or gcc
      --cxx CXX             C++ compiler to use, default is CXX environment variable or g++
      --fc FC               Fortran compiler to use, defualt is FC environment variable or gfortran
      --download-dir DOWNLOAD_DIR
                            Download location of tarballs
      --extract-dir EXTRACT_DIR
                            Extract and Build location
      --build-shared        Build shared libraries (Default)
      --build-static        Build static libraries
      -j JOBS, --jobs JOBS  Number of parallel jobs
      --enable-parmetis     Build ParMETIS library, (Default, check license requirements)
      --disable-parmetis    Disable ParMETIS library
      --cmake-dir CMAKE_DIR
                            System location of package cmake
      --openmpi-dir OPENMPI_DIR
                            System location of package openmpi
      --hdf5-dir HDF5_DIR   System location of package hdf5
      --pnetcdf-dir PNETCDF_DIR
                            System location of package pnetcdf
      --netcdf-dir NETCDF_DIR
                            System location of package netcdf
      --fmt-dir FMT_DIR     System location of package fmt
      --seacas-dir SEACAS_DIR
                            System location of package seacas
      --bison-dir BISON_DIR
                            System location of package bison
      --flex-dir FLEX_DIR   System location of package flex
      --openblas-dir OPENBLAS_DIR
                            System location of package openblas
      --metis-dir METIS_DIR
                            System location of package metis
      --parmetis-dir PARMETIS_DIR
                            System location of package parmetis
      --scotch-dir SCOTCH_DIR
                            System location of package scotch
      --arpack_ng-dir ARPACK_NG_DIR
                            System location of package arpack-ng
      --scalapack-dir SCALAPACK_DIR
                            System location of package scalapack
      --mumps-dir MUMPS_DIR
                            System location of package mumps
      --superlu_dist-dir SUPERLU_DIST_DIR
                            System location of package superlu_dist
      --suitesparse-dir SUITESPARSE_DIR
                            System location of package suitesparse
      --trilinos-dir TRILINOS_DIR
                            System location of package trilinos
      --petsc-dir PETSC_DIR
                            System location of package petsc
      --petsc_complex-dir PETSC_COMPLEX_DIR
                            System location of package petsc-complex
      --omega_h-dir OMEGA_H_DIR
                            System location of package omega-h
      --sparse-dir SPARSE_DIR
                            System location of package sparse