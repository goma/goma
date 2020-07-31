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


#include "wr_soln.h"
#define _WR_SOLN_C
#include "goma.h"

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void write_solution(char output_file[],             /* name EXODUS II file */
                    double resid_vector[],          /* Residual vector */
                    double x[],                     /* soln vector */
                    double **x_sens_p,              /* sensitivity vectors */
                    double x_old[],                 /* soln vector at previous time step */
                    double xdot[],                  /* dx/dt */
                    double xdot_old[],              /* dx/dt at previous time step */
                    int tev,                        /* total elem variables
                                                       (normal) */
                    int tev_post,                   /* additional post process
                                                       elem vars */
                    double *gv,                     /* global variable values */
                    struct Results_Description *rd, /* for post process vars
                                                       (rf_io_structs.h) */
                    int *gindex,
                    int *p_gsize,
                    double *gvec,
                    double ***gvec_elem, /* array vars [eb_index][ev_index][elem] */
                    int *nprint,         /* counter for time step number*/
                    dbl delta_t,         /* time step size */
                    dbl theta,           /* time integration parameter */
                    dbl time_value,      /* current time value */
                    dbl *x_pp,           /* post proc vars for export */
                    Exo_DB *exo,
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
  int i, i_post, step = 0;
#ifdef DEBUG
  static char *yo = "write_solution";
  fprintf(stderr, "%s: begins\n", yo);
#endif

  /* First nodal quantities */
  for (i = 0; i < rd->TotalNVSolnOutput; i++) {
    extract_nodal_vec(x, rd->nvtype[i], rd->nvkind[i], rd->nvmatID[i], gvec, exo, FALSE,
                      time_value);
    step = (*nprint) + 1;
    wr_nodal_result_exo(exo, output_file, gvec, i + 1, step, time_value);
  }

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

    post_process_nodal(x, x_sens_p, x_old, xdot, xdot_old, resid_vector, step, &time_value, delta_t,
                       theta, x_pp, exo, dpi, rd, output_file);
#ifdef DEBUG
    fprintf(stderr, "%s: done w/ post_process_nodal\n", yo);
#endif

    /*
     *  Write out time derivatives if requested
     */
    if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
      for (i = 0; i < rd->TotalNVSolnOutput; i++) {
        i_post = rd->TotalNVSolnOutput + rd->TotalNVPostOutput + i;
        extract_nodal_vec(xdot, rd->nvtype[i_post], rd->nvkind[i_post], rd->nvmatID[i], gvec, exo,
                          TRUE, time_value);
        wr_nodal_result_exo(exo, output_file, gvec, i_post + 1, *nprint + 1, time_value);
      }
    }
  }

  /* Now element quantities */
  for (i = 0; i < tev; i++) {
    bool is_P1 = FALSE;
    int dof = 0;
    for (int mn = 0; mn < upd->Num_Mat; mn++) {
      if (pd_glob[mn]->i[rd->evtype[i]] == I_P1) {
        dof = MAX(getdofs(type2shape(exo->eb_elem_itype[mn]), I_P1), dof);
        is_P1 = TRUE;
      }
    }
    if (is_P1) {
      for (int k = 0; k < dof; k++) {
        extract_elem_vec(x, i, rd->evtype[i], gvec_elem, exo, k);
        step = (*nprint) + 1;
        wr_elem_result_exo(exo, output_file, gvec_elem, i, step, time_value, rd);
        i++;
      }
    } else {
      extract_elem_vec(x, i, rd->evtype[i], gvec_elem, exo, 0);
      step = (*nprint) + 1;
      wr_elem_result_exo(exo, output_file, gvec_elem, i, step, time_value, rd);
    }
  }
  /* Finally, global values */

  wr_global_result_exo(exo, output_file, step, rd->ngv, gv);

#ifdef DEBUG
  fprintf(stderr, "%s: done with regular element vars; start tev_post\n", yo);
#endif

  /* Add additional user-specified post processing variables */
  if (tev_post > 0) {
    step = (*nprint) + 1;
#ifdef DEBUG
    fprintf(stderr, "%s: start post_process_elem\n", yo);
#endif

    post_process_elem(x, x_old, xdot, xdot_old, resid_vector, tev, tev_post, gvec_elem, step,
                      &time_value, delta_t, exo, dpi, rd);
#ifdef DEBUG
    fprintf(stderr, "%s: done w/ post_process_elem\n", yo);
#endif

    /* Write out time derivatives if requested */
    if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
      for (i = 0; i < tev; i++) {
        i_post = tev_post + i;
        extract_elem_vec(xdot, i_post, rd->evtype[i_post], gvec_elem, exo, 0);
        wr_elem_result_exo(exo, output_file, gvec_elem, i_post, *nprint + 1, time_value, rd);
      }
    }
  }
}
/*****************************************************************************/
/*  END of file wr_soln.c  */
/*****************************************************************************/
