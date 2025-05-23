/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: rf_mp.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

/*
 * Rename to avoid Aztec name conflicts.
 */

#ifndef GOMA_RF_MP_H
#define GOMA_RF_MP_H

#include <mpi.h>

/*
 * Parameters describing the characteristics of the Parallel Machine
 * and utility vectors used during MPI communications to monitor status
 * of sends and recvs (particularly, nonblocking)
 *
 * These variables are declared here for general purpose usage. They
 * are defined in main.c, where the arrays are also allocated and freed.
 */

extern int Num_Proc; /* Total number of processors */

extern int ProcID; /* This processor's number */

extern int Dim; /* Dimension of logical hypercube */

extern MPI_Request *Request;

extern MPI_Status *Status;

extern int Num_Requests;

#endif /* end of GOMA_RF_MP_H */
