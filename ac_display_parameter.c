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
 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "std.h"

#include "exo_struct.h"
#include "dpi.h"

#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "rf_solver_const.h"
#include "rf_solver.h"

#include "rf_masks.h"

#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp.h"
#include "mm_mp_structs.h"

#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_eh.h"

#include "exo_struct.h"		/* defn of Exo_DB */
#include "dp_types.h"

#include "sl_util.h"		/* defines sl_init() */

extern void display_parameterC
PROTO((double, 			/* PARAMETER VALUE */
       double*, 		/* UNKNOWN VECTOR */
       double*, 		/* UNKNOWN_DOT VECTOR */
       double, 			/* STEP */
       Comm_Ex *,		/* cx  - array of communications structures */
       Exo_DB *,		/* exo - ptr to finite element mesh database */
       Dpi *));			/* dpi - ptr to distributed processing info */

void
display_parameterC(double lambda, /* PARAMETER VALUE */
		   double *x, 	 /* UNKNOWN VECTOR */
		   double *xdot,  /* UNKNOWN_DOT VECTOR */
		   double delta_s,/* STEP */
		   Comm_Ex *cx,	 /* array of communications structures */
		   Exo_DB *exo,	 /* ptr to the finite element mesh database */
		   Dpi *dpi)	 /* distributed processing information */
{
  int ic;
  int mn;
  int ibc, idf;
  
#ifdef DEBUG
  static char yo[]="display_parameterC";
#endif

  /*
   * 		BEGIN EXECUTION
   */

#ifdef DEBUG
  fprintf(stderr, "display_parameter() begins...\n");
#endif

  if (cont->upType == 1) {

    ibc = cont->upBCID;
    idf = cont->upDFID;

    printf("\n New BC[%4d] DF[%4d] = %10.6e\n", ibc, idf, lambda);

  }
  else
  if (cont->upType == 2) {

    mn = cont->upMTID;
    ic = cont->upMDID;

    printf("\n New MT[%4d] MP[%4d] = %10.6e\n", mn+1, cont->upMPID, lambda);
 
  }

}/* END of routine display_parameterC  */
/*******************************************************************************/



