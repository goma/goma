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
 

#include "wr_soln.h"

#include <stdio.h>

#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_more_utils.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io_structs.h"
#include "wr_exo.h"
#include "dpi.h"
#include "exo_struct.h"
#include "std.h"

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void 
write_solution(char output_file[], double resid_vector[], double x[], double **x_sens_p, double x_old[], double xdot[],
               double xdot_old[], int tev, int tev_post, double *gv, struct Results_Description *rd, double *gvec,
               double ***gvec_elem, int *nprint, dbl delta_t, dbl theta, dbl time_value, dbl *x_pp, Exo_DB *exo,
               Dpi *dpi)
    /*************************************************************************
     *
     * write_solution():
     *
     * This routine gathers both the solution and the post-processed
     * variables and then writes them one at a time to the output file.
     * The description of what it writes is taken from the
     * Results_Description structure in the parameter list.
     *
     *************************************************************************/
{
  int i, i_post, step=0;
#ifdef DEBUG
  static char *yo="write_solution";
  fprintf(stderr, "%s: begins\n", yo);
#endif

  /* First nodal quantities */
  for (i = 0; i < rd->TotalNVSolnOutput; i++) {
    extract_nodal_vec(x, rd->nvtype[i], rd->nvkind[i], rd->nvmatID[i],
		      gvec, exo, FALSE, time_value);
    step = (*nprint)+1;
    wr_nodal_result_exo(exo, output_file, gvec, i+1, step, time_value);
  }

  /* Special case for global post processing, usually file output */
#ifdef DEBUG
  fprintf(stderr, "%s: done with regular nodal vars; start global\n", yo);
#endif

  /* Special case for global post processing, special case file output for now*/
  post_process_global(x, exo, dpi, time_value);

#ifdef DEBUG
  fprintf(stderr, "%s: done with global; start tnv_post\n", yo);
#endif
  /*
   *  Add additional user-specified post processing variables
   */
  if (rd->TotalNVPostOutput > 0) {
    step = (*nprint) + 1;
#ifdef DEBUG
    fprintf(stderr, "%s: start post_process_nodal\n", yo);
#endif

    post_process_nodal(x, x_sens_p, x_old, xdot, xdot_old, resid_vector, 
		       step, &time_value, delta_t, theta, x_pp,
		       exo, dpi, rd, output_file, 0);
#ifdef DEBUG
    fprintf(stderr, "%s: done w/ post_process_nodal\n", yo);
#endif

    /*
     *  Write out time derivatives if requested
     */
    if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
      for (i = 0; i < rd->TotalNVSolnOutput; i++) {
	i_post = rd->TotalNVSolnOutput + rd->TotalNVPostOutput + i;
	extract_nodal_vec(xdot, rd->nvtype[i_post], rd->nvkind[i_post],
			  rd->nvmatID[i], gvec, exo, TRUE, time_value);
	wr_nodal_result_exo(exo, output_file, gvec, i_post + 1, 
			    *nprint+1, time_value);
      }
    }
  }


  /* Now element quantities */
  for(i = 0; i < tev; i++) {
    extract_elem_vec(x, i, rd->evtype[i], gvec_elem, exo);
    step = (*nprint)+1;
    wr_elem_result_exo(exo, output_file, gvec_elem, i, step, 
		       time_value, rd);
  }
  /* Finally, global values */
  
  wr_global_result_exo( exo, output_file, step, rd->ngv, gv );
	
#ifdef DEBUG
  fprintf(stderr, "%s: done with regular element vars; start tev_post\n", yo);
#endif

  /* Add additional user-specified post processing variables */
  if (tev_post > 0) {
      step = (*nprint) + 1;
#ifdef DEBUG
      fprintf(stderr, "%s: start post_process_elem\n", yo);
#endif

      post_process_elem(x, x_old, xdot, xdot_old, resid_vector, tev, tev_post, 
			gvec_elem, step, &time_value, delta_t, exo, dpi, rd);
#ifdef DEBUG
      fprintf(stderr, "%s: done w/ post_process_elem\n", yo);
#endif

      /* Write out time derivatives if requested */
      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
	for (i = 0; i < tev; i++) {
	  i_post = tev_post + i;
	  extract_elem_vec(xdot, i_post, rd->evtype[i_post], gvec_elem, exo);
	  wr_elem_result_exo(exo, output_file, gvec_elem, i_post,
			     *nprint+1, time_value, rd);
	}
      }
  }
}


void
write_solution_segregated(char output_file[], double **resid_vector, double **x, double **x_old, double **xdot,
                          double **xdot_old, int *tev_post, double *gv, struct Results_Description **rd,
                          double **gvec, int *nprint, dbl delta_t, dbl theta, dbl time_value, dbl *x_pp,
                          Exo_DB *exo, Dpi *dpi)
{
  int i, step=0;
  int i_post;
#ifdef DEBUG
  static char *yo="write_solution";
  fprintf(stderr, "%s: begins\n", yo);
#endif

  /* First nodal quantities */
  int offset = 0;
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    /* Special case for global post processing, usually file output */
    post_process_global(x[pg->imtrx], exo, dpi, time_value);

    if (pg->imtrx > 0) {
      offset += rd[pg->imtrx -1]->TotalNVSolnOutput + rd[pg->imtrx -1]->TotalNVPostOutput;
      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
	offset += rd[pg->imtrx -1]->TotalNVSolnOutput;
      }
    }
    for (i = 0; i < rd[pg->imtrx]->TotalNVSolnOutput; i++) {
      extract_nodal_vec(x[pg->imtrx], rd[pg->imtrx]->nvtype[i], rd[pg->imtrx]->nvkind[i],
          rd[pg->imtrx]->nvmatID[i], gvec[pg->imtrx], exo, FALSE, time_value);
      step = (*nprint) + 1;
      wr_nodal_result_exo(exo, output_file, gvec[pg->imtrx], offset + i + 1, step, time_value);
    }
  }

#ifdef DEBUG
  fprintf(stderr, "%s: done with regular nodal vars; start tnv_post\n", yo);
#endif

  /*
   *  Add additional user-specified post processing variables
   */
  offset = 0;
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    if (pg->imtrx > 0) {
      offset += rd[pg->imtrx -1]->TotalNVSolnOutput + rd[pg->imtrx -1]->TotalNVPostOutput;
      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
	offset += rd[pg->imtrx -1]->TotalNVSolnOutput;
      }
    }
    if (rd[pg->imtrx]->TotalNVPostOutput > 0) {
      step = (*nprint) + 1;
  #ifdef DEBUG
      fprintf(stderr, "%s: start post_process_nodal\n", yo);
  #endif

      post_process_nodal(x[pg->imtrx], NULL, x_old[pg->imtrx], xdot[pg->imtrx],
			 xdot_old[pg->imtrx], resid_vector[pg->imtrx], step, &time_value,
			 delta_t, theta, x_pp, exo, dpi, rd[pg->imtrx], output_file, offset);
  #ifdef DEBUG
      fprintf(stderr, "%s: done w/ post_process_nodal\n", yo);
  #endif

      /*
       *  Write out time derivatives if requested
       */

    }
  }

  offset = 0;
  if (tev_post != NULL) {
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

      step = (*nprint) + 1;

      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
	if (pg->imtrx > 0) {
	  offset += rd[pg->imtrx-1]->TotalNVSolnOutput + rd[pg->imtrx-1]->TotalNVPostOutput + rd[pg->imtrx-1]->TotalNVSolnOutput;
	}
	for (i = 0; i < rd[pg->imtrx]->TotalNVSolnOutput; i++) {
	  i_post = rd[pg->imtrx]->TotalNVSolnOutput + rd[pg->imtrx]->TotalNVPostOutput + i;
	  extract_nodal_vec(xdot[pg->imtrx], rd[pg->imtrx]->nvtype[i_post], rd[pg->imtrx]->nvkind[i_post],
			    rd[pg->imtrx]->nvmatID[i], gvec[pg->imtrx], exo, TRUE, time_value);
	  wr_nodal_result_exo(exo, output_file, gvec[pg->imtrx], offset + i_post + 1,
			      *nprint+1, time_value);
	}
      }
    }
  }
//
  /* Now element quantities */
//  for(i = 0; i < tev; i++) {
//    extract_elem_vec(x, i, rd->evtype[i], gvec_elem, exo);
//    step = (*nprint)+1;
//    wr_elem_result_exo(exo, output_file, gvec_elem, i, step,
//                       time_value, rd);
//  }
  /* Finally, global values */

  wr_global_result_exo( exo, output_file, step, rd[0]->ngv, gv );

#ifdef DEBUG
  fprintf(stderr, "%s: done with regular element vars; start tev_post\n", yo);
#endif

  /* Add additional user-specified post processing variables */
//  if (tev_post > 0) {
//      step = (*nprint) + 1;
//#ifdef DEBUG
//      fprintf(stderr, "%s: start post_process_elem\n", yo);
//#endif
//
//      post_process_elem(x, x_old, xdot, xdot_old, resid_vector, tev, tev_post,
//                        gvec_elem, step, &time_value, delta_t, exo, dpi, rd);
//#ifdef DEBUG
//      fprintf(stderr, "%s: done w/ post_process_elem\n", yo);
//#endif
//
//      /* Write out time derivatives if requested */
//      if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
//        for (i = 0; i < tev; i++) {
//          i_post = tev_post + i;
//          extract_elem_vec(xdot, i_post, rd->evtype[i_post], gvec_elem, exo);
//          wr_elem_result_exo(exo, output_file, gvec_elem, i_post,
//                             *nprint+1, time_value, rd);
//        }
//      }
//  }
}
/*****************************************************************************/
/*  END of file wr_soln.c  */
/*****************************************************************************/
