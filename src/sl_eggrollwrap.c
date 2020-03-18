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
 
/*
 * $Id: sl_eggrollwrap.c,v 5.2 2007-09-18 18:53:47 prschun Exp $
 */

#ifdef USE_RCSID
static const char rcs_id[] = "$Id: sl_eggrollwrap.c,v 5.2 2007-09-18 18:53:47 prschun Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ac_stability.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "mm_more_utils.h"
#include "mm_post_proc.h"
#include "rd_exo.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_umf.h"
#include "std.h"
#include "wr_dpi.h"
#include "wr_exo.h"

/* Routines that handle the eigensolver.
 *
 * Linear stability analysis
 * Solve J z = t M z
 *
 * where 
 *
 * input:
 * J = jacobian matrix 
 * M = mass or overlap matrix
 *
 * output:
 * z = eigenvectors
 * t = eigenvalues
 *
 * Friendly warning:
 * Do not edit this unless you know what you are doing!!
 *
 * Originally written by Ian Gates
 * pre-CVS modification history:
 *   - Sep 24, 1997, first checkin
 *   - Feb 98 -> Oct 98, another checkin
 *   - Jan 13, 2000, MMH rearranged and conformed to Goma style.
 */
void
eggrollwrap(int *istuff,	/* info for eigenvalue extraction */
	    dbl *dstuff,	/* info for eigenvalue extraction */
	    
	    int *ija,		/* column pointer array */
	    dbl *jac,		/* nonzero array */
	    dbl *mas,		/* nonzero array - same structure 
				   as jac[] (ija[]) */
	    
	    dbl *x,		/* Value of the solution vector */
	    char *ExoFileOut,	/* Name of exoII output file */
	    int prob_type,
	    dbl delta_t,	/* time step size */
	    dbl theta,		/* variable time integration parameter
				   explicit (theta = 1) to 
				   implicit (theta = 0) */
	    dbl *x_old,		/* Value of the old solution vector */
	    dbl *xdot,		/* Value of xdot predicted for new 
				   solution */
	    dbl *xdot_old,      /* dx/dt at previous time step */
	    dbl *resid_vector,
	    int *converged,	/* whether the Newton has converged */
	    int *nprint,	/* counter for time step number */
	    int tnv,		/* number of nodal results */
	    int tnv_post,	/* number of post processing results */
	    struct Results_Description *rd,
	    int *gindex,
	    int *p_gsize,
	    dbl	*gvec,
	    dbl time_value,
	    Exo_DB *exo,	/* ptr to finite element mesh db */
	    int Num_Proc,	/* number of processors used */
	    Dpi	*dpi)		/* ptr to distributed processing info */
{
  int 
    i, j, ic, 
    nj, nnz_j, 
    first_linear_solver_call, Factor_Flag, matr_form, 
    error, rcflag, action, 
    ev_n, ev_jac, 
    filter, 
    mm, max_itr,
    nev_want, nev_found, lead, 
    /*    read_form, soln_tech, push_mode, */
    push_mode, 
    init_shft, recycle;
  dbl
    stol, ivector,
    dwork[20]; 
  dbl 
    *ev_e, *ev_i, *ev_r, *ev_x,  
    *v1, *v2, 
    *mat, 
    **evect, **schur;
  char save_ExoFileOut[MAX_FNL];

  static int UMF_system_id;	/* Used to uniquely identify the
				 * explicit fill system to solve from
				 * the other UMF systems. */
  /* Initialize
   */
  ic = error = rcflag = action = 0;
  ev_jac = 0;
  matr_form = 1;

  /* Set values
   */
  mm         = istuff[0];
  nj         = istuff[1];
  nnz_j      = istuff[2];
  filter     = istuff[3];
  recycle    = istuff[4];
  nev_want   = istuff[6];
  init_shft  = istuff[7];
  max_itr    = istuff[8];
  push_mode  = istuff[9];
  stol       = dstuff[0];
  ivector    = dstuff[3];

  printf(" Initializing variables and allocating space... ");

  /* Allocate spectrum storage
   */
  ev_e = Dvector_birth(mm+5);
  ev_i = Dvector_birth(mm+5);
  ev_r = Dvector_birth(mm+5);
  ev_x = Dvector_birth(mm+5);

  /* Set initial (real) shifts
   */
  ev_n = init_shft;
  vcopy(init_shft, &ev_r[0], 1.0, &dstuff[10]);

  /* Allocate auxiliary work vectors
   */
  mat = Dvector_birth(nnz_j+5);

  /* Allocate eigenvectors and schur storage
   */
  i = nj+5;
  j = mm+5; 
  evect = Dmatrix_birth(j, i);
  schur = Dmatrix_birth(j, i);

  /* Allocate reverse communication vectors
   */
  v1 = Dvector_birth(nj+5);
  v2 = Dvector_birth(nj+5);

  /* Check for something that seems to make no difference if it's on,
   * except for occasionally causing seg faults... */
  if(recycle != 0)
    EH(-1, "Eigen recycle currently doesn't work, turn it off.");

  /* Set initial vector
   */
  vinit(nj, &v1[0], 0.5);

  /* GEVP solution
   */
  ic = 0;
  first_linear_solver_call = +1;
  do {
    ic++;
    /* printf("ic = %d\n", ic); fflush(stdout); */
    gevp_solver_rc(nj, mm, max_itr, stol, filter, &ev_n, &ev_r[0],
		   &ev_i[0], &ev_e[0], &ev_x[0], &lead, ev_jac,
		   nev_want, &nev_found, schur, evect, recycle,
		   ivector, 0, &rcflag, &action, &dwork[0], &v1[0],
		   &v2[0]);
    /* printf("action = %d\n", action); fflush(stdout); */
    switch (action)
      {
      case  0: /* All done */
	break;  
      case  1: /* v2 = J*v1 */
	MV_MSR(&nj, &ija[0], &jac[0], &v1[0], &v2[0]);
	break;  
      case  2: /* v2 = M*v1 */
	MV_MSR(&nj, &ija[0], &mas[0], &v1[0], &v2[0]);
	break;  
      case  3: /* inv(J-sM) */
	/* Shift matrix step */
	v2sum(nnz_j, &mat[0], 1.0, &jac[0], -dwork[0], &mas[0]);

	/* Invert step - get LU for later */
	if(first_linear_solver_call == 1)
	  {
	    Factor_Flag = -2;
	    UMF_system_id = -1;
	  }
	else
	  Factor_Flag = -1;

	/*
	printf("Calling SL_UMF, first_linear_solver_call = %d, Factor_Flag = %d\n",
	       first_linear_solver_call, Factor_Flag); fflush(stdout); 
	*/

	UMF_system_id = SL_UMF(UMF_system_id,
			       &first_linear_solver_call, 
			       &Factor_Flag, 
			       &matr_form, 
			       &nj, 
			       &nnz_j, 
			       &ija[0], 
			       &ija[0], 
			       &mat[0], 
			       &v1[0], 
			       &v2[0]);
	first_linear_solver_call = 0;
	break;
      case  4: /* v2 = inv(J-sM)*M*v1 */
	Factor_Flag = 3;
	if(first_linear_solver_call)
	  EH(-1, "Tried to transform eigenvectors before a solve!");
	gevp_transformation(UMF_system_id, first_linear_solver_call,
			    Factor_Flag, matr_form, 1, nj, nnz_j,
			    &ija[0], &jac[0], &mas[0], &mat[0],
			    /*			  soln_tech, */
			    &v2[0], &v1[0], dwork[0], dwork[1]);
	break;
      default:
	EH(-1, "Uh-oh!  I shouldn't be here!");
	break;
      } /* switch(action) */
    if (ic > 10000)
      error = 1;
  } while ((rcflag != 0) && (error == 0));

  /* Error check
   */
  if (error == 1)
    {
      puts(" E: Too many iterations.  Escape.  ");
      exit(-1);
    }

  /* De-allocate solver storage
   */  

  first_linear_solver_call = -1;

  /* MMH sez: If first_linear_solver_call == -1, then we want to
   * deallocate memory so we shouldn't be trying to solve anything!
   * This was FMR'ing b/c SL_UMF was being called with
   * first_linear_solver_call = -1, and Factor_Flag = 3.  Bad.
   */
  Factor_Flag = -3;

  UMF_system_id = SL_UMF(UMF_system_id,
			 &first_linear_solver_call, 
			 &Factor_Flag, 
			 &matr_form, 
			 &nj, 
			 &nnz_j, 
			 ija, 
			 ija, 
			 mat, 
			 &v1[0], 
			 &v2[0]);  
  
  /* Display results 
   */
  printf("\n-------------------------------------------------------------------------------\n");
  if(Linear_Stability == LSA_3D_OF_2D)
    printf("NORMAL MODE WAVE NUMBER = %g\n", LSA_3D_of_2D_wave_number);
  printf(" Eigensolver required %d iterations.\n",ic);
  printf(" Found %d converged eigenvalues.\n", nev_found);
  printf(" Leading Eigenvalue  = % 10.6e%+10.6e i RES = % 10.6e\n", 
	 ev_r[lead], ev_i[lead], ev_e[lead]);
  printf("    Real           Imag           RES\n");
  for (i=0;i<nev_found;i++)
    printf(" % 10.6e %+10.6e i % 10.6e\n", ev_r[i], ev_i[i], ev_e[i]);

  /* MMH: I know this is stupid, but the filename for the "regular"
   * Exodus output is a global variable!!!  It is required in
   * post_process_nodal().  I swap it out here, and will swap it back
   * when we're done with LSA.  Why don't I just overwrite it
   * completely you may ask?  Well, I don't know if and/or when the
   * code will continue to do something useful after LSA.  If it ever
   * does, then it would probably like to know what the correct output
   * filename is.  Kinda like camping: Leave with what you came in
   * with.  */
  strncpy(save_ExoFileOut, ExoFileOut, MAX_FNL-1);

  /* Write results to file (exoII format)
   */
  printf(" push_mode                          = %12d  \n", push_mode);
  if (push_mode > 0)
    {
      puts(" Writing modes to file ...");
      /* Write to exo file
       * Each mode is written as a "time step" solution into exoII DB
       */
      for(i = 0; i < push_mode; i++)
	{
	  printf("\t\t Mode %4d ...", i);
	  if(LSA_3D_of_2D_wave_number == -1.0)
	    sprintf(ExoFileOut, "LSA_%d_of_%d_%s", i + 1, push_mode,
		    save_ExoFileOut);
	  else
	    sprintf(ExoFileOut, "LSA_%d_of_%d_wn=%g_%s", i + 1, push_mode,
		    LSA_3D_of_2D_wave_number, save_ExoFileOut);

	  /* Replicate basic mesh info */
	  one_base(exo);
	  wr_mesh_exo(exo, ExoFileOut, 0);
	  zero_base(exo);
	  wr_result_prelim_exo(rd, exo, ExoFileOut, NULL);
	  /* Update exo file for distributed problem info 
	   */
	  if (Num_Proc > 1) {
	    wr_dpi(dpi, ExoFileOut, 0);
	  }
	  for (j = 0; j < tnv; j++) {
	    extract_nodal_vec(&evect[i][0], rd->nvtype[j], rd->nvkind[j], 
			      rd->nvmatID[j], gvec, exo, FALSE, time_value);
	    wr_nodal_result_exo(exo, ExoFileOut, gvec, j+1, 1, 
				time_value);
	  }

	  /*
	   *  Add additional user-specified post processing variables 
	   */
	  if (tnv_post > 0) {
	    post_process_nodal(&evect[i][0], NULL, x_old, xdot, xdot_old,
			       resid_vector, 1, &time_value, delta_t, 0.0,
                               NULL, exo, dpi, rd, ExoFileOut, 0);
	  }
	  printf(" recorded.\n");
	}
    }
  /* MMH: See comments above. */
  strncpy(ExoFileOut, save_ExoFileOut, MAX_FNL);

  /* De-allocate work vectors
   */
  printf("Deallocating memory ... ");
  i = nj+5;
  j = mm+5;
  Dmatrix_death(schur, j, i);
  Dmatrix_death(evect, j, i);
  Dvector_death(&v2[0], nj+5);
  Dvector_death(&v1[0], nj+5);
  Dvector_death(&mat[0], nnz_j+5);
  Dvector_death(&ev_e[0], mm+5);
  Dvector_death(&ev_i[0], mm+5);
  Dvector_death(&ev_r[0], mm+5);
  Dvector_death(&ev_x[0], mm+5);
  printf("done.\n");
}
