/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 

#undef linux
#undef sun
#undef dec_osf1
#undef _AIX
#undef sgi
#undef hpux

/* DEBUG mode(s) */
#undef DEBUG

/* MPI installation */
#undef PARALLEL
#undef SERIAL
#undef HAVE_MPI
#undef MPI

/* Frontal solver */
#undef HAVE_FRONT

/* Sparse library */
#undef HAVE_SPARSE

/* UMFPACK */
#undef HAVE_UMFPACK

/* BLAS */
#undef HAVE_Y12M

/* BLAS */
#undef GOMA_HAVE_BLAS

/* LAPACK */
#undef GOMA_HAVE_LAPACK

/* Aztec */
#undef HAVE_AZTEC

/* ARPACK */
#undef HAVE_ARPACK
#ifdef HAVE_ARPACK
#define EIGEN_SERIAL
#endif

/* PARPACK */
#undef HAVE_PARPACK
#ifdef HAVE_PARPACK
#define EIGEN_PARALLEL
#undef EIGEN_SERIAL
#endif

