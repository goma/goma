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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_particles.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_update_parameter.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_rotate.h"
#include "bc_special.h"
#include "bc_surfacedomain.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "std.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_out_hkm.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h"
#endif


/***************************************************************************/
/************* GLOBAL FUNCTIONS IN THIS FILE     ***************************/
/***************************************************************************
*
*  user_out ()           void                  rf_solve, rf_transient, etc.
*
****************************************************************************/

static void PRNT_SOLN_VALS(char *soln_name, double l_max, double g_max,
                           int i_max, double l_min, double g_min,
                           int i_min, double g_avg);

static int GET_PRINT_PROC_NUM(const double, const double);
static void RETN_COORD(int node, double *x, double *y, double *z );

static void SOLN_VALUES(int, int, double [], double *,
                        double *, int *, double *, double *,
                        int *, double *, const int);

static void SOLN_TO_MOLES(double [], double [], double [], double [],
			  MATRL_PROP_STRUCT *);
static void GET_SOLN_LASTSPECIES(double[], double [], MATRL_PROP_STRUCT *);

#ifdef USE_CHEMKIN
#ifdef SENKIN_OUTPUT
static void WRITE_SENKIN_FILE(int, int, int,  double [],  double,
                              int,  MATRL_PROP_STRUCT *matID_prop  );
#endif
#endif

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void 
usr_out_hkm(int status, double time, double dt, double *soln)     
   
  /*************************************************************************
  * user_out_hkm:
  *
  * User-defined function that allows for the user to do any kind of output
  * desired.  It is called at the times the user probably wants to do output.
  * The situation for the call is represented by the status variable.
  *
  * Values for status variable:
  *-------------------------------
  * < 0   = Some sort of error condition has occurred.
  *   0   = Initial conditions
  *   1   = Final conditions, i.e., a successful run has completed
  *   2   = A successful intermediate time step has occurred.
  *   3   = A successful intermediate step in the non-linear solver
  *         has occurred.
  *
  * The current solution and time is passed in to the routine.  The current
  * time step number is also passed into the routine.  This refers to the
  * exodus time step number, i.e., the initial conditions are called time
  * step number 1.
  **************************************************************************/
{
  int mn;        /* Material ID Number */
  static int    time_step_num = 0;

  double          *soln_mole = NULL; /* solution vector with mole fractions
                                        instead of mass fractions          */
  double          *soln_lastSpecies = NULL;/* Temporary vector to store 
					      the calculated value of the
					      mf of last species */
  double    *soln_mole_lastSpecies = NULL;
  double           Tmin, Tmax, g_Tmin, g_Tmax, Tavg;
  int              k, num_owned_unknowns;
  int              iTmax, iTmin, var;
  char             char_str[72];
  struct Problem_Description *pd_local;
  MATRL_PROP_STRUCT *matID_prop;
  /*  double  depo_rate, non_uniformity; */
  int Debug_Flag = 3;

#ifdef USE_CHEMKIN
#ifdef SENKIN_OUTPUT
  int              SenkinNode = 0;
#endif
#endif
  
  time_step_num++;

  /************ Print out a Banner at the start of User I/O **************/
  /*
   * Limit when output is to occur - output for final, every successful time
   * step, and for initial conditions
   */
  
  if (Debug_Flag <= 2) return;

  /* Number of unknowns owned by this processor -- rf_fem.h */
  num_owned_unknowns  = NumUnknowns[pg->imtrx];
  
#ifdef PARALLEL
  (void) MPI_Barrier(MPI_COMM_WORLD);
#endif  
  if (ProcID == 0) {
    printf("\n\n\t");
    fprint_line(stdout,"u", 70);
    printf("\tPrintout of the user defined solution statistics\n\t");
    printf("\tTime = %g, time step num = %d\n\t", time, time_step_num);
    fprint_line(stdout,"u", 70);
  }
  
  /*
   * Loop over the number of materials
   *      -> We will print out totals for each material
   */
  for (mn = 0; mn < upd->Num_Mat; mn++) {
    pd_local = pd_glob[mn];
    matID_prop = mp_glob[mn];
      
    /************** Fluid Mechanics Printouts  *****************************/
    var = VELOCITY1;
    if (pd_local->v[0][var]) {
      SOLN_VALUES(var, 0, soln, &Tmax, &g_Tmax, &iTmax, &Tmin,
		  &g_Tmin, &iTmin, &Tavg, mn);

      PRNT_SOLN_VALS("X Velocity", Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		     Tavg);
    }
    var = VELOCITY2;
    if (pd_local->v[0][var]) {
      SOLN_VALUES(var, 0, soln, &Tmax, &g_Tmax, &iTmax, &Tmin,
		  &g_Tmin, &iTmin, &Tavg, mn);

      PRNT_SOLN_VALS("Y Velocity", Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		     Tavg);
    }
    var = VELOCITY3;
    if (pd_local->v[0][var]) {
      SOLN_VALUES(var, 0, soln, &Tmax, &g_Tmax, &iTmax, &Tmin,
		  &g_Tmin, &iTmin, &Tavg, mn);

      PRNT_SOLN_VALS("Z Velocity", Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		     Tavg);
    }
    var = PRESSURE;
    if (pd_local->v[0][var]) {
      SOLN_VALUES(var, 0, soln, &Tmax, &g_Tmax, &iTmax, &Tmin,
		  &g_Tmin, &iTmin, &Tavg, mn);

      PRNT_SOLN_VALS("Pressure", Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		     Tavg);
    }

    /************ Heat Transfer Printout ***********************************/

    var = TEMPERATURE;
    if (pd_local->v[0][var]) {
      SOLN_VALUES(var, 0, soln, &Tmax, &g_Tmax, &iTmax, &Tmin,
		  &g_Tmin, &iTmin, &Tavg, mn);

      PRNT_SOLN_VALS("Temperature", Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		     Tavg);
    }

    /*********** Mass Transfer Temporary Printouts *************************/

    var = MASS_FRACTION;
    if (pd_local->v[0][var]) {
      /*
       *   Write out a SENKIN file for the first node on the 0th processor,
       *   if chemkin was used as the default database.
       */
#ifdef USE_CHEMKIN
#ifdef SENKIN_OUTPUT
      if (ProcID == 0) {
	if (matID_prop->DefaultDatabase == DB_CHEMKIN_MAT) {
	  WRITE_SENKIN_FILE(status, 1,  SenkinNode, soln, time,
			    time_step_num, matID_prop);
	}
      }
#endif
#endif
         
      /*
       * For CHEMKIN material types, translate the unknowns to mole fractions
       * before printing them out.  We will want to do this for multi-species
       * material types also, for the case where the mass fractions sum to
       * one. However, we can't be sure that this is the case.  Therefore, 
       * the conversion to mole fractions is indeterminable.
       */

	
      for (k = 0; k < matID_prop->Num_Species_Eqn; k++) {
	SOLN_VALUES(MASS_FRACTION, k, soln, &Tmax, &g_Tmax,
		    &iTmax, &Tmin, &g_Tmin, &iTmin, &Tavg, mn);
	(void) sprintf(char_str, "Mass Fraction for Sp# %2d (%.16s)",k+1,
		       matID_prop->Species_Names[k]);
	PRNT_SOLN_VALS(char_str, Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
		       Tavg);
      }

      /*
       * Calculate the mass fraction of the last species in the
       * mechanism. And then print it out, like all of the rest.
       * Malloc a temporary vector, the size of the total solution vector
       * on the processor. We will store the last solution for the
       * the last species in the first mass fraction slot in
       * this temporary vector -> wasteful but what the heck.
       */
      if (matID_prop->Num_Species_Eqn < matID_prop->Num_Species) {
	soln_lastSpecies  = alloc_dbl_1(num_owned_unknowns, 0.0);
	GET_SOLN_LASTSPECIES(soln, soln_lastSpecies, matID_prop);
	SOLN_VALUES(MASS_FRACTION, 0, soln_lastSpecies, &Tmax, &g_Tmax,
		    &iTmax, &Tmin, &g_Tmin, &iTmin, &Tavg, mn);
	sprintf(char_str, "Mass Fraction for Sp# %2d (%.16s)", k+1,
		matID_prop->Species_Names[k]);
	PRNT_SOLN_VALS(char_str, Tmax, g_Tmax, iTmax, Tmin,
		       g_Tmin, iTmin, Tavg);
      }
#ifdef PARALLEL
#endif
      /*
       * Translate the solution from mass to mole fractions. Then,
       * print it out.
       */
      if (pd_local->Species_Var_Type == SPECIES_MASS_FRACTION) {
	/*
	 * Create another solution vector which is mole fraction based.
	 * Don't bother to include external nodes.
	 */
	soln_mole             = alloc_dbl_1(num_owned_unknowns, 0.0);
	soln_mole_lastSpecies = alloc_dbl_1(num_owned_unknowns, 0.0);
	SOLN_TO_MOLES(soln, soln_mole, soln_lastSpecies,
		      soln_mole_lastSpecies, matID_prop);
#ifdef PARALLEL
#endif
	for (k = 0; k < matID_prop->Num_Species_Eqn; k++) {
	  SOLN_VALUES(MASS_FRACTION, k, soln_mole, &Tmax, &g_Tmax,
		      &iTmax, &Tmin, &g_Tmin, &iTmin, &Tavg, mn);
	  sprintf(char_str, "Mole Fraction for Sp# %2d (%.16s)", k+1,
		  matID_prop->Species_Names[k]);
	  PRNT_SOLN_VALS(char_str, Tmax, g_Tmax, iTmax, Tmin, g_Tmin, iTmin,
			 Tavg);
	}
	if (matID_prop->Num_Species_Eqn != matID_prop->Num_Species) {
	  SOLN_VALUES(MASS_FRACTION, 0, soln_mole_lastSpecies, &Tmax, &g_Tmax,
		      &iTmax, &Tmin, &g_Tmin, &iTmin, &Tavg, mn);
	  k = matID_prop->Num_Species_Eqn;
	  sprintf(char_str, "Mole Fraction for Sp# %2d (%.16s)", k+1,
		  matID_prop->Species_Names[k]);
	  PRNT_SOLN_VALS(char_str, Tmax, g_Tmax, iTmax, Tmin,
			 g_Tmin, iTmin, Tavg);
	}
      }
#ifdef PARALLEL
#endif
      /*
     * Free Memory
     */
      safer_free((void **) &soln_mole);
      safer_free((void **) &soln_lastSpecies);
      safer_free((void **) &soln_mole_lastSpecies);
    }
    /************* Print a Banner Line at the End of Output ********************/
    if (ProcID == 0) {
      printf("\n\t"); fprint_line(stdout,"u",50);
      printf("\tEnd of of the user defined solution statistics\n\t");
      fprint_line(stdout,"u",50);
    }
  }
} /************ END of user_out () *******************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void PRNT_SOLN_VALS(char *soln_name, double l_max, double g_max,
                           int i_max, double l_min, double g_min,
                           int i_min, double g_avg)
{
  int    print_proc;
  double Xcoor, Ycoor, Zcoor;

  /*
   * Write a header out, and write average value out Process the information
   * for the maximum value
   */

  fflush(stdout);
  print_proc = GET_PRINT_PROC_NUM(l_max, g_max);
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (ProcID == print_proc) {
    printf("\t %s:\n", soln_name);
    printf("\t  Average value = %g\n", g_avg);

    RETN_COORD(i_max, &Xcoor, &Ycoor, &Zcoor);
    printf("\t  Max value = %g\tat Xcoor: %g\tYcoor: %g", g_max, Xcoor, Ycoor);
    if (Num_Dim > 2)
	printf("\tZcoor: %g", Zcoor);
    printf("\n");
    fflush(stdout);
  }
  
  /*
   * Process the information for the minimum value
   */
  print_proc = GET_PRINT_PROC_NUM(l_max, g_max);
  if (ProcID == print_proc) {  
    RETN_COORD(i_min, &Xcoor, &Ycoor, &Zcoor);
    printf("\t  Min value = %g\tat Xcoor: %g\tYcoor: %g", g_min,
	   Xcoor, Ycoor);
    if (Num_Dim > 2)
	printf("\tZcoor: %g", Zcoor);
    printf("\n");
    fflush(stdout);
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
} /* PRNT_SOLN_VALS */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int GET_PRINT_PROC_NUM(const double num, const double g_num)

{
  int proc_temp, print_proc;
  if (num == g_num)   proc_temp = ProcID;
  else                proc_temp = Num_Proc;
  print_proc = gmin_int(proc_temp);
  return print_proc;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void SOLN_VALUES(int var_type, int sub_index, double soln[],
			double *l_max, double *g_max,  int *i_max,
                        double *l_min,  double *g_min,  int *i_min,
			double *g_avg, const int mn)
    
    /************************************************************************
     *
     * SOLN_VALUES()
     *
     * Find the maximum and minimum values of a component of the soln
     * vector. Also, finds the average value. mn is the material number.
     * If it is negative, then the average is taken over all materials.
     ************************************************************************/
{
  double  s_max = -DBL_MAX, s_min = DBL_MAX, s_sum = 0.0, s, g_sum;
  int     ebi, i_eqn, num_nodes_mn = 0, g_solncomp_nodes, index, n, i;
  int     e_start, e_end, ielem, ielem_type, num_local_nodes, ndof = 0;
  int *i_been_there;

  int num_owned_nodes = DPI_ptr->num_internal_nodes +
                        DPI_ptr->num_boundary_nodes;
  

  *i_max = *i_min = 0;

  /*
   *  Set the number of nodes equal to the number of owned nodes on this
   *  processor
   */
  i_been_there = alloc_int_1(num_owned_nodes, 0);

  for (ebi = 0; ebi < EXO_ptr->num_elem_blocks; ebi++) {
    if (Matilda[ebi] == mn) {
      e_start = EXO_ptr->eb_ptr[ebi];
      e_end   = EXO_ptr->eb_ptr[ebi+1];
      for (ielem = e_start; ielem < e_end; ielem++) {
	ielem_type = Elem_Type(EXO_ptr, ielem);
	num_local_nodes = elem_info(NNODES, ielem_type);
	index      = Proc_Connect_Ptr[ielem];
	for (n = 0; n < num_local_nodes; n++) {
	  i = Proc_Elem_Connect[index++];
	  if (Nodes[i]->Type.Owned) {
	    if (! i_been_there[i]) {
	      i_been_there[i] = 1;
	
	      /*
	       * Find the max, min, and sum of the soln component on the
	       * current processor
	       * - Note: not all soln components are on all processors,
	       *   so there are additional complications
	       */
              ndof = 0;
	      if ((i_eqn =
		   Index_Solution(i, var_type, sub_index, ndof, mn, pg->imtrx)) >= 0) {
		num_nodes_mn++;
		s      = soln[i_eqn];
		s_sum += fabs(s);

		if (s > s_max) {
		  s_max  = s;
		  *i_max = i;
		}
		if (s < s_min) {
		  s_min  = s;
		  *i_min = i;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /*
   * Set the local processor maximum and minimum values
   */
  *l_max = s_max;
  *l_min = s_min;

  /* Compute the global max and min values */

  if (Num_Proc > 1) {
#ifdef PARALLEL
    MPI_Allreduce(&s_max, g_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&s_min, g_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    /*
     * Find the total number of nodes containing the soln component
     */
    MPI_Allreduce(&num_nodes_mn, &g_solncomp_nodes, 1, MPI_INT,
			 MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&s_sum, &g_sum, 1, MPI_DOUBLE,
			 MPI_SUM, MPI_COMM_WORLD);
#endif
  } else {
    *g_max = s_max;
    *g_min = s_min;
    g_solncomp_nodes = num_nodes_mn;
    g_sum = s_sum;
  }
 
  /*
	     * Compute the average value of the Component
	     */
  if (g_solncomp_nodes > 0) {
    *g_avg = g_sum / g_solncomp_nodes;
  } else {
    *g_avg = 0.0;
  }
  safer_free((void **) &i_been_there);
} /* SOLN_VALUES */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void RETN_COORD(int node, double *x, double *y, double *z)

{
  *x = Coor[0][node];
  *y = Coor[1][node];
  if (Num_Dim > 2)  *z = Coor[2][node];
  else              *z = 0.0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void SOLN_TO_MOLES(double soln[], double soln_mole[],
			  double soln_lastSpecies[],
			  double soln_mole_lastSpecies[],
			  MATRL_PROP_STRUCT *matID_prop)

    /************************************************************************
     *
     * SOLN_TO_MOLES()
     *
     * Translate the solution from mass to mole fractions
     ************************************************************************/
{
  int    i_eqn, j_eqn, k, var, mn;
  double sumwt;
  int  node, ndof = 0;
  int num_owned_nodes = (DPI_ptr->num_internal_nodes +
			 DPI_ptr->num_boundary_nodes);
  mn = matID_prop->MatID;
  for (node = 0; node < num_owned_nodes; node++) {
    if ((i_eqn =
	 Index_Solution(node, MASS_FRACTION, 0, ndof, mn, pg->imtrx)) >= 0) { 
      for (k = 0, sumwt = 0.0; k < matID_prop->Num_Species_Eqn; k++) {
	sumwt +=  soln[i_eqn + k] / matID_prop->molecular_weight[k];
      }
      if (matID_prop->Num_Species >  matID_prop->Num_Species_Eqn) {
	k = matID_prop->Num_Species_Eqn;
	sumwt += soln_lastSpecies[i_eqn] / matID_prop->molecular_weight[k];
      }
      for (k = 0; k < matID_prop->Num_Species_Eqn; k++) {
	soln_mole[i_eqn + k] = soln[i_eqn + k] /
            (sumwt * matID_prop->molecular_weight[k]);
      }
      if (matID_prop->Num_Species > matID_prop->Num_Species_Eqn) {
	k = matID_prop->Num_Species_Eqn;
	soln_mole_lastSpecies[i_eqn] = soln_lastSpecies[i_eqn]/
                     (sumwt * matID_prop->molecular_weight[k]);
      }
    }
  }

  /*
   * Take advantage of the extra space and check that mole fractions sum to one
   */

  var = TEMPERATURE;
  if (pd_glob[mn]->v[pg->imtrx][var]) {
    for (node = 0; node < num_owned_nodes; node++) {
      if ((i_eqn =
	   Index_Solution(node, MASS_FRACTION, 0, ndof, mn, pg->imtrx)) >= 0) {
        if ((j_eqn = Index_Solution(node, TEMPERATURE, 0, ndof, mn, pg->imtrx)) >= 0) {
          for (k = 0, sumwt = 0.0; k < matID_prop->Num_Species_Eqn; k++) {
            sumwt += soln_mole[i_eqn + k];
          }
	  if (matID_prop->Num_Species >  matID_prop->Num_Species_Eqn) {
	     sumwt += soln_mole_lastSpecies[i_eqn];
	  }
          soln_mole[j_eqn] = sumwt;
        }
      }
    }
  }
} /* SOLN_TO_MOLES */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void GET_SOLN_LASTSPECIES(double soln[], double soln_lastSpecies[],
                                 MATRL_PROP_STRUCT *matID_prop)
  
    /***********************************************************************
     *
     *  GET_SOLN_LASTSPECIES:
     *
     *     Get the solution for the last species in the mechanism.
     ***********************************************************************/
{
  int i_eqn, k;
  int node, ndof = 0;
  int num_owned_nodes = (DPI_ptr->num_internal_nodes +
			 DPI_ptr->num_boundary_nodes);
  int mn = matID_prop->MatID;
  ndof = mn;
  if (soln_lastSpecies == NULL) return;
  if (matID_prop->Num_Species_Eqn < matID_prop->Num_Species) {
    for (node = 0; node < num_owned_nodes; node++) {
      if ((i_eqn =
	   Index_Solution(node, MASS_FRACTION, 0, ndof, mn, pg->imtrx)) >= 0) {
	soln_lastSpecies[i_eqn] = 1.0;
	for (k = 0; k < matID_prop->Num_Species_Eqn; k++) {
	  soln_lastSpecies[i_eqn] -= soln[i_eqn + k];
	}
      }
    }
  } else {
    for (node = 0; node < num_owned_nodes; node++) {
      if ((i_eqn =
	   Index_Solution(node, MASS_FRACTION, 0, ndof, mn, pg->imtrx)) >= 0) {
	soln_lastSpecies[i_eqn] = soln[i_eqn];
      }
    }
  }
} /********************** GET_SOLN_LASTSPECIES *****************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
#ifdef USE_CHEMKIN
#ifdef SENKIN_OUTPUT

static void WRITE_SENKIN_FILE(int    status,
                              int    mesh_id,
                              int node,
                              double soln[],
                              double time,
                              int    time_step_num,
                              MATRL_PROP_STRUCT *matID_prop
                              )

    /***********************************************************************
     *
     * WRITE_SENKIN_FILE():
     *
     * This routine writes a SENKIN binary output file, which describes the 
     * time behavior of the solution at a particular node, node.
     ***********************************************************************/
{
  static int    LUNIT;
  static int    LOUT = 6;
  static int    LSENS;
  static int    NSYS, KK, II, ICASE, NEQ;
  static double TSTOP, TLIM, TOLS[4], VOL;
  int           indx_MF, indx_T, var, k;
  int           ndof = 0;
  int           mn = matID_prop->MatID;
  double        T, pressureAmbient, sum;
  double       *MF_vector = NULL;
  
  extern FSUB_TYPE write_header_(int *, int *, int *, int *, int *, int *,
                                 int *, int *, double *, double *, double *,
                                 int *, double *, CK_NAME *, int );
  extern FSUB_TYPE write_record_(int *, int *, double *, double *, int *,
			         int *, int *, int *, double *,
                                 double *, double *, double *);

  if (status == 0) {

    /*
     * Open the file for writing using a fortran function
     */

    LUNIT = senkin_open_();
    if (LUNIT == 0) {
      (void) fprintf(stderr,"save.dat failed to open\n");
      exit (-1);
    }

    /*
     * Output of the Header information
     */

    LSENS = 0;
    KK    = matID_prop->Num_Species;
    NSYS  = 1 + KK;
    II    = Ck.IiGas;

    /*
     * Use the constant pressure, calculate the temperature case for senkin,
     * even if the temperature is not a variable
     */

    ICASE = 1;
    NEQ   = KK + 1;
    TSTOP = TimeMax;
    TLIM  = 400.;
    TOLS[0] = 1.0E-3;
    TOLS[1] = 1.0E-8;
    TOLS[2] = TOLS[0];
    TOLS[3] = TOLS[1];

    /*
     * Call a FORTRAN function that opens the file, save.dat, and writes out
     * the necessary header information
     */

    (void) write_header_(&LOUT, &LUNIT, &LSENS, &NSYS, &KK, &II, &ICASE, 
                         &NEQ, &TSTOP, &TLIM, TOLS, Ck.Ickwrk_p, Ck.Rckwrk_p,
                         Ck.Cckwrk_p, 16);
  }

  /*
   * Write a valid record for status = 1 or 2
   */

  VOL = 1.0;

  if ((indx_MF = Index_Solution(node, MASS_FRACTION, 0, ndof, mn, pg->imtrx)) < 0) {
    (void) fprintf(stderr, "error indx_MF is bad\n");
    exit(-1);
  }

  var = TEMPERATURE;
  if (pd_glob[0]->v[pg->imtrx][var]) {
    if ((indx_T =
	 Index_Solution(mesh_id, node, TEMPERATURE, 0, mn, pg->imtrx)) >= 0) {
      T = soln[indx_T];
    }
#ifdef DEBUG
    else {
      EH(-1,"error indx_T is bad\n");
    }
#endif
  }
  else {
    T = matID_prop->reference[TEMPERATURE];
  }

  pressureAmbient = 1.0133E6;
  /*
   * Call the fortran subroutine that writes out a record of data
   */
  MF_vector = (double *) array_alloc(1, KK, sizeof(double));

  sum = 0.0;
  for (k = 0; k < matID_prop->Num_Species_Eqn; k++) {
    MF_vector[k] = soln[indx_MF + k];
    sum += soln[indx_MF + k];
  }
  if (matID_prop->Num_Species_Eqn < matID_prop->Num_Species) {
         MF_vector[matID_prop->Num_Species_Eqn] = 1.0 - sum; 
  }

  (void) write_record_(&LUNIT, &ICASE, &time, MF_vector, &NSYS,
		       &LSENS, &II, Ck.Ickwrk_p, Ck.Rckwrk_p, &T,
                       &pressureAmbient, &VOL);

  safer_free((void **) &MF_vector);

  /*
   * Close the output file, if we are at the end of the calculation
   */

  if (status == 1) {
    file_close_(&LUNIT);
  }

} /* WRITE_SENKIN_FILE */

#endif
#endif
/***************************************************************************/

