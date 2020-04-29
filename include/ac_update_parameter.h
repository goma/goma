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
 
#ifndef GOMA_AC_UPDATE_PARAMETER_H
#define GOMA_AC_UPDATE_PARAMETER_H

#include "dp_types.h"
#include "exo_struct.h"
#include "dpi.h"
#include "ac_hunt.h"

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_UPDATE_PARAMETER_C
#define EXTERN
#
#endif

#ifndef GOMA_AC_UPDATE_PARAMETER_C
#define EXTERN extern
#endif

EXTERN void update_parameterC
(int,                     /* CONDITION NUMBER */
       double, 			/* PARAMETER VALUE */
       double*, 		/* UNKNOWN VECTOR */
       double*, 		/* UNKNOWN_DOT VECTOR */
       double*, 		/* x_AC VECTOR */
       double, 			/* STEP */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_parameterTP
(int,                     /* CONDITION NUMBER */
       double, 			/* PARAMETER VALUE */
       double*, 		/* UNKNOWN VECTOR */
       double*, 		/* UNKNOWN_DOT VECTOR */
       double*, 		/* x_AC VECTOR */
       double, 			/* STEP */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_parameterHC
(int,                     /* HC ID */
       double, 			/* PARAMETER VALUE */
       double*, 		/* UNKNOWN VECTOR */
       double*, 		/* UNKNOWN_DOT VECTOR */
       double*, 		/* x_AC VECTOR */
       double, 			/* STEP */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_parameterS
( double, 		/* PARAMETER VALUE */
        double* , 		/* UNKNOWN VECTOR */
     	double* ,		/* UNKNOWN_DOT VECTOR */
	  int ,			/*  type of sensitivity variable(BC or MT) */
	  int ,			/*  variable id BCID or MT# */
	  int ,			/* data float id or matl prop id */
	  int ,			/* data float id or matl prop id */
	  Comm_Ex *,		/* array of communications structures */
	  Exo_DB *,		/* ptr to the finite element mesh database */
	  Dpi *);		/* distributed processing information */

EXTERN void update_BC_parameter
(double, 			/* PARAMETER VALUE */
       int,     		/* BCID */
       int,     		/* DFID */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_AC_parameter
(double, 			/* PARAMETER VALUE */
       int,     		/* BCID */
       int,     		/* DFID */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_MT_parameter
(double, 			/* PARAMETER VALUE */
       int, 			/* Material number index */
       int,     		/* Material property tag */
       int,     		/* Material property tag subindex */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void update_UM_parameter
(double,                  /* PARAMETER VALUE */
       int,                     /* Material number index */
       int,     		/* Material property tag */
       int,     		/* Material property tag subindex */
       Comm_Ex *,               /* cx  - array of communications structures */
       Exo_DB *,                /* exo - ptr to finite element mesh database */
       Dpi *);                 /* dpi - ptr to distributed processing info */

EXTERN void retrieve_parameterS
(	double* , /* PARAMETER VALUE */
        double* , 	 /* UNKNOWN VECTOR */
     	double* ,  /* UNKNOWN_DOT VECTOR */
	  int ,/*  type of sensitivity variable(BC or MT) */
	  int ,		/*  variable id BCID or MT# */
	  int ,		/* data float id or matl prop id */
	  int ,		/* matl prop id for UM */
	  Comm_Ex *,	 /* array of communications structures */
	  Exo_DB *,	 /* ptr to the finite element mesh database */
	  Dpi *);	 /* distributed processing information */

EXTERN void retrieve_BC_parameter
(double*, 			/* PARAMETER VALUE */
       int,     		/* BCID */
       int,     		/* DFID */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void retrieve_MT_parameter
(double*, 			/* PARAMETER VALUE */
       int, 			/* Material number index */
       int,     		/* Material property tag */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void retrieve_AC_parameter
(double*, 			/* PARAMETER VALUE */
       int,     		/* BCID */
       int,     		/* DFID */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

EXTERN void retrieve_UM_parameter
(double*,                  /* PARAMETER VALUE */
       int,                     /* Material number index */
       int,     		/* Material property tag */
       int,     		/* Material property tag subindex */
       Comm_Ex *,               /* cx  - array of communications structures */
       Exo_DB *,                /* exo - ptr to finite element mesh database */
       Dpi *);                 /* dpi - ptr to distributed processing info */

#endif /* GOMA_AC_UPDATE_PARAMETER_H */
