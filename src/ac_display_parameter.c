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
 

#include <stdio.h>

#include "exo_struct.h"
#include "dpi.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "dp_types.h"

extern void display_parameterC
(double, 			/* PARAMETER VALUE */
       double*, 		/* UNKNOWN VECTOR */
       double*, 		/* UNKNOWN_DOT VECTOR */
       double, 			/* STEP */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *);			/* dpi - ptr to distributed processing info */

void
display_parameterC(double lambda, /* PARAMETER VALUE */
		   double *x, 	 /* UNKNOWN VECTOR */
		   double *xdot,  /* UNKNOWN_DOT VECTOR */
		   double delta_s,/* STEP */
		   Comm_Ex *cx,	 /* array of communications structures */
		   Exo_DB *exo,	 /* ptr to the finite element mesh database */
		   Dpi *dpi)	 /* distributed processing information */
{
  int mn;
  int ibc, idf;
  

  /*
   * 		BEGIN EXECUTION
   */


  if (cont->upType == 1) {

    ibc = cont->upBCID;
    idf = cont->upDFID;

    printf("\n New BC[%4d] DF[%4d] = %10.6e\n", ibc, idf, lambda);

  }
  else
  if (cont->upType == 2) {

    mn = cont->upMTID;

    printf("\n New MT[%4d] MP[%4d] = %10.6e\n", mn+1, cont->upMPID, lambda);
 
  }

}/* END of routine display_parameterC  */
/*******************************************************************************/



