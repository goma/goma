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

/* dp_vif.c -- distributed processing utilities for virtual input file
 *
 * [1]	Load up the initialization information structure that embodies the
 *      essential information from the GOMA ASCII input file and material
 *      files. That is, if we're on the main I/O processor.
 *
 * [2]  Broadcast it from this proc. Receive on all other procs.
 *
 * [3]  Once received, unload the structure into the standard extern variables
 *      and other places.
 *
 * [4]  Three routines and three stages. Raven, Ark, Dove. This accomodates
 *      various dynamically allocated data.
 *
 * [5]  Include every file that declares some important global external variable
 *      that needs to be referenced...
 *
 *
 * Created: 1997/05/19 11:03 MDT pasacki@sandia.gov
 *
 * Revised: 1997/06/19 08:48 MDT pasacki@sandia.gov
 */

#include <stdio.h>
#include <stdlib.h>

#include "ac_particles.h"
#include "ac_stability_util.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_chemkin.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mpi.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_mp.h"
#include "rf_solver.h"
#include "std.h"

/* #include "mm_names.h" -- use extern defs from rf_bc_const.h for these vars */

#include "rf_bc.h"

#define GOMA_DP_VIF_C

#ifndef CDIM
#define CDIM 3
#endif

/*
 * Setup types for stuff that is of known size and setup stuff used to size
 * other stuff...
 */

extern DDD Noahs_Dove;

extern DDD Noahs_Raven;

extern DDD Noahs_Ark;

/*
 * Noah's Raven is the first look. Here, set up information that is needed
 * so that there is a proper landing place for most of the data. Also,
 * this helps economize the data transfer somewhat for some problems where
 * meaningful data only partially fills big static data structures.
 */

void noahs_raven(void) {
#ifdef PARALLEL
  DDD n;

  /*
   *  Also send the upd structure information at this time
   *  -> Note, every item in this structure is of fixed length
   *     Therefore, just broadcast the structure.
   */

  ddd_add_member2(upd, 1, sizeof(UPD_STRUCT));
  ddd_set_commit2();

  Noahs_Raven = ddd_alloc();
  n = Noahs_Raven;

  ddd_add_member(n, &Decompose_Flag, 1, MPI_INT);
  ddd_add_member(n, &Decompose_Type, 1, MPI_INT);
  ddd_add_member(n, &Skip_Fix, 1, MPI_INT);

  ddd_add_member(n, &Num_Var_Init, 1, MPI_INT);
  ddd_add_member(n, &Num_Var_Bound, 1, MPI_INT);
  ddd_add_member(n, &Num_Var_External, 1, MPI_INT);
  ddd_add_member(n, &TimeIntegration, 1, MPI_INT);
  ddd_add_member(n, &Num_BC, 1, MPI_INT);
  ddd_add_member(n, &num_new_BC_Desc, 1, MPI_INT);
  ddd_add_member(n, &Num_ROT, 1, MPI_INT);
  ddd_add_member(n, &Num_Interpolations, 1, MPI_INT);
  ddd_add_member(n, &CoordinateSystem, 1, MPI_INT);
  ddd_add_member(n, &len_u_post_proc, 1, MPI_INT);
  ddd_add_member(n, &num_AC_Tables, 1, MPI_INT);
  ddd_add_member(n, &num_BC_Tables, 1, MPI_INT);
  ddd_add_member(n, &num_MP_Tables, 1, MPI_INT);
  ddd_add_member(n, &num_ext_Tables, 1, MPI_INT);
  ddd_add_member(n, &Continuation, 1, MPI_INT);
  ddd_add_member(n, &nAC, 1, MPI_INT);
  ddd_add_member(n, &nCC, 1, MPI_INT);
  ddd_add_member(n, &nTC, 1, MPI_INT);
  ddd_add_member(n, &nHC, 1, MPI_INT);
  ddd_add_member(n, &nUC, 1, MPI_INT);
  ddd_add_member(n, &nUTC, 1, MPI_INT);
  ddd_add_member(n, &nEQM, 1, MPI_INT);
  ddd_add_member(n, &Linear_Stability, 1, MPI_INT);
  ddd_add_member(n, &LSA_number_wave_numbers, 1, MPI_INT);
  ddd_add_member(n, &nn_post_fluxes, 1, MPI_INT);
  ddd_add_member(n, &nn_post_fluxes_sens, 1, MPI_INT);
  ddd_add_member(n, &nn_error_metrics, 1, MPI_INT);
  ddd_add_member(n, &nn_post_data, 1, MPI_INT);
  ddd_add_member(n, &nn_post_data_sens, 1, MPI_INT);
  ddd_add_member(n, &nn_volume, 1, MPI_INT);
  ddd_add_member(n, &nn_global, 1, MPI_INT);
  ddd_add_member(n, &nn_average, 1, MPI_INT);
  ddd_add_member(n, &Chemkin_Needed, 1, MPI_INT);
  ddd_add_member(n, &efv->ev, 1, MPI_INT);
  ddd_add_member(n, &LOCA_UMF_ID, 1, MPI_INT);
  ddd_add_member(n, &Use_Level_Set, 1, MPI_INT);
  ddd_add_member(n, &Use_Phase_Field, 1, MPI_INT);
  ddd_add_member(n, &Use_DG, 1, MPI_INT);
  ddd_add_member(n, &Do_Overlap, 1, MPI_INT);
  ddd_add_member(n, &Particle_Dynamics, 1, MPI_INT);
  ddd_add_member(n, &Particle_Number_Sample_Types, 1, MPI_INT);
  ddd_add_member(n, &Particle_Number_PBCs, 1, MPI_INT);
  ddd_add_member(n, &Num_Var_LS_Init, 1, MPI_INT);
  ddd_add_member(n, &TFMP_GAS_VELO, 1, MPI_INT);
  ddd_add_member(n, &TFMP_LIQ_VELO, 1, MPI_INT);
  ddd_add_member(n, &TFMP_INV_PECLET, 1, MPI_INT);
  ddd_add_member(n, &TFMP_KRG, 1, MPI_INT);
  ddd_add_member(n, &GomaPetscOptionsStrLen, 1, MPI_INT);

  ddd_set_commit(n);

#endif
  return;
}

/*
 * raven_post_alloc() -- allocate space after 1st round of communication
 *
 * Once the raven has carried this essential information to other processors,
 * allocate some space so the ark can land.
 *
 * Created: 1997/06/19 09:02 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void raven_landing(void) {
  int i;
  int m;

  if (ProcID == 0)
    return;

  /*
   * These are allocations that processor 0 does during the course of
   * reading the goma input file. Since the other processors do not read
   * the input file, this allocation is done here separately, once the
   * required size information has been communicated via the raven.
   */

  upd->turbulent_info = calloc(1, sizeof(turbulent_information));

  /*
   * Instead of communicating all efv, just this bit, the remainder in
   * the ark.
   */

  efv->Num_external_field = Num_Var_External;

  /*
   * Gratuitous replication of this is achieved on other processors here.
   */

  for (i = 0; i < MAX_NUMBER_MATLS; i++) {
    pd_glob[i]->TimeIntegration = TimeIntegration;
  }

  if (Num_BC > 0) {
    BC_Types =
        (struct Boundary_Condition *)array_alloc(1, Num_BC, sizeof(struct Boundary_Condition));
  }

  /* Broadcast input wave number array for 3D of 2D LSA */
  if (LSA_number_wave_numbers > 0 &&
      (Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE)) {
    LSA_wave_numbers = (double *)array_alloc(1, LSA_number_wave_numbers, sizeof(double));
  } else {
    LSA_wave_numbers = NULL;
  }

  if (GomaPetscOptionsStrLen > 0) {
    GomaPetscOptions = calloc(GomaPetscOptionsStrLen, sizeof(char));
  }

  /*
   * Create landing pads.
   */

  /*
   * Create landing pad for new BC_Description structures
   */

  if (num_new_BC_Desc > 0) {
    new_BC_Desc = (struct BC_descriptions **)array_alloc(1, num_new_BC_Desc,
                                                         sizeof(struct BC_descriptions *));
    for (i = 0; i < num_new_BC_Desc; i++) {
      new_BC_Desc[i] = (struct BC_descriptions *)calloc(1, sizeof(struct BC_descriptions));
    }
  }
  /*
   * Create landing pads for TABLE_AC Data_Table structures
   */

  for (i = 0; i < num_AC_Tables; i++) {
    AC_Tables[i] = smalloc(sizeof(struct Data_Table));
  }

  /*
   * Create landing pads for TABLE_BC Data_Table structures
   */

  for (i = 0; i < num_BC_Tables; i++) {
    BC_Tables[i] = smalloc(sizeof(struct Data_Table));
  }

  /*
   * Create landing pads for TABLE_MP Data_Table structures
   */

  for (i = 0; i < num_MP_Tables; i++) {
    MP_Tables[i] = smalloc(sizeof(struct Data_Table));
  }

  /*
   * Create landing pads for external tables Data_Table structures
   */

  for (i = 0; i < num_ext_Tables; i++) {
    ext_Tables[i] = smalloc(sizeof(struct Data_Table));
  }

  /*
   * Rotation condition preliminaries...
   */

  if (Num_ROT > 0) {
    ROT_Types = (struct Rotation_Specs *)calloc(Num_ROT, sizeof(struct Rotation_Specs));
  }

  /*
   * Material property preliminaries that are already known...
   */

  for (m = 0; m < upd->Num_Mat; m++) {
    pd_glob[m]->CoordinateSystem = CoordinateSystem;
  }

  /*
   * User defined post processing with variable number of user provided
   * constants.
   */

  if (len_u_post_proc > 0) {
    u_post_proc = (dbl *)array_alloc(1, len_u_post_proc, sizeof(dbl));
  }

  /*
   * flux post processing with variable number of user provided
   * constants.
   */

  if (nn_post_fluxes > 0) {
    pp_fluxes = (struct Post_Processing_Fluxes **)array_alloc(
        1, nn_post_fluxes, sizeof(struct Post_Processing_Fluxes *));
    for (i = 0; i < nn_post_fluxes; i++) {
      pp_fluxes[i] =
          (struct Post_Processing_Fluxes *)array_alloc(1, 1, sizeof(struct Post_Processing_Fluxes));
    }
  }

  if (nn_post_fluxes_sens > 0) {
    pp_fluxes_sens = (struct Post_Processing_Fluxes_Sens **)array_alloc(
        1, nn_post_fluxes_sens, sizeof(struct Post_Processing_Fluxes_Sens *));

    for (i = 0; i < nn_post_fluxes_sens; i++) {
      pp_fluxes_sens[i] = (struct Post_Processing_Fluxes_Sens *)array_alloc(
          1, 1, sizeof(struct Post_Processing_Fluxes_Sens));
    }
  }

  if (nn_post_data > 0) {
    pp_data = (struct Post_Processing_Data **)array_alloc(1, nn_post_data,
                                                          sizeof(struct Post_Processing_Data *));
    for (i = 0; i < nn_post_data; i++) {
      pp_data[i] =
          (struct Post_Processing_Data *)array_alloc(1, 1, sizeof(struct Post_Processing_Data));
    }
  }

  if (nn_post_data_sens > 0) {
    pp_data_sens = (struct Post_Processing_Data_Sens **)array_alloc(
        1, nn_post_data_sens, sizeof(struct Post_Processing_Data_Sens *));

    for (i = 0; i < nn_post_data_sens; i++) {
      pp_data_sens[i] = (struct Post_Processing_Data_Sens *)array_alloc(
          1, 1, sizeof(struct Post_Processing_Data_Sens));
    }
  }
  /*
   * flux post processing sensitivities  with variable number of user provided
   * constants.
   */
  if (nn_volume > 0) {
    pp_volume = (struct Post_Processing_Volumetric **)array_alloc(
        1, nn_volume, sizeof(struct Post_Processing_Volumetric *));

    for (i = 0; i < nn_volume; i++) {
      pp_volume[i] = (struct Post_Processing_Volumetric *)array_alloc(
          1, 1, sizeof(struct Post_Processing_Volumetric));
    }
  }

  if (nn_global > 0) {
    pp_global = (struct Post_Processing_Global **)array_alloc(
        1, nn_global, sizeof(struct Post_Processing_Global *));
    for (i = 0; i < nn_global; i++) {
      pp_global[i] =
          (struct Post_Processing_Global *)array_alloc(1, 1, sizeof(struct Post_Processing_Global));
    }
  }

  if (nn_average > 0) {
    pp_average = (struct Post_Processing_Averages **)array_alloc(
        1, nn_average, sizeof(struct Post_Processing_Averages *));
    for (i = 0; i < nn_average; i++) {
      pp_average[i] = (struct Post_Processing_Averages *)array_alloc(
          1, 1, sizeof(struct Post_Processing_Averages));
    }
  }
  /*
   * Zienkewicz-Zhu error measures.
   */
  /*
    if ( nn_error_metrics > 0 )
      {
        pp_error_data = (struct Post_Processing_Error **)
         array_alloc(1, nn_error_metrics, sizeof(struct Post_Processing_Error *));
        for ( i = 0; i < nn_error_metrics; i++)
          {
            pp_error_data[i] = (struct Post_Processing_Error *)
                      array_alloc(1, 1, sizeof(struct Post_Processing_Error));
          }

      }
  */

  /*
   * Augmenting condition preliminaries...
   */

  if (nAC > 0) {
    augc = alloc_struct_1(struct AC_Information, nAC);
  } else {
    augc = NULL;
  }

  /*
   * Continuation condition preliminaries...
   */

  if (nCC > 0) {
    cpcc = alloc_struct_1(struct Continuation_Conditions, nCC);
  } else {
    cpcc = NULL;
  }

  /*
   * Turning point continuation condition preliminaries...
   */

  if (nTC > 0) {
    tpcc = (struct Continuation_Conditions *)array_alloc(1, nTC,
                                                         sizeof(struct Continuation_Conditions));
  } else {
    tpcc = NULL;
  }

  /*
   * User continuation condition preliminaries...
   */

  if (nUC > 0) {
    cpuc =
        (struct User_Continuation_Info *)array_alloc(1, nUC, sizeof(struct User_Continuation_Info));
  } else {
    cpuc = NULL;
  }

  /*
   * User turning point continuation condition preliminaries...
   */

  if (nUTC > 0) {
    tpuc = (struct User_Continuation_Info *)array_alloc(1, nUTC,
                                                        sizeof(struct User_Continuation_Info));
  } else {
    tpuc = NULL;
  }

  /*
   * Hunting condition preliminaries...
   */

  if (nHC > 0) {
    hunt = (struct HC_Information *)array_alloc(1, nHC, sizeof(struct HC_Information));
  } else {
    hunt = NULL;
  }

  /*
   * Make the Level Set Structure
   **/
  if (Use_Level_Set) {
    ls = (struct Level_Set_Data *)calloc(1, sizeof(struct Level_Set_Data));
    ls->embedded_bc = NULL;
    ls->init_surf_list = NULL;
    lsi = (struct Level_Set_Interface *)calloc(1, sizeof(struct Level_Set_Interface));
    zero_lsi();
    zero_lsi_derivs();

    for (m = 0; m < upd->Num_Mat; m++) {
      mp_glob[m]->mp2nd = (SECOND_LS_PHASE_PROP_STRUCT *)alloc_void_struct_1(
          sizeof(SECOND_LS_PHASE_PROP_STRUCT), 1);
    }

  } else {
    ls = NULL;
    lsi = NULL;
  }

  /*
   * Make the Phase field Structure
   **/
  if (Use_Phase_Field) {
    pfd = alloc_struct_1(struct Phase_Function_Data, 1);
    pfd->num_phase_funcs = Use_Phase_Field;
    pfd->ls = (struct Level_Set_Data **)alloc_ptr_1(pfd->num_phase_funcs);

    for (i = 0; i < pfd->num_phase_funcs; i++) {
      pfd->ls[i] = alloc_struct_1(struct Level_Set_Data, 1);
      pfd->ls[i]->embedded_bc = NULL;
      pfd->ls[i]->init_surf_list = NULL;
    }

  } else {
    pfd = NULL;
  }

  if (Particle_Dynamics) {
    if (Particle_Number_Sample_Types) {
      Particle_Number_Samples_Existing =
          (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Number_Samples = (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Number_Output_Variables =
          (int *)array_alloc(1, Particle_Number_Sample_Types, sizeof(int));
      Particle_Output_Variables = (particle_variable_s **)array_alloc(
          1, Particle_Number_Sample_Types, sizeof(particle_variable_s *));
      Particle_Filename_Template = (particle_filename_s *)array_alloc(
          1, Particle_Number_Sample_Types, sizeof(particle_filename_s));
    } else {
      Particle_Number_Samples_Existing = NULL;
      Particle_Number_Samples = NULL;
      Particle_Number_Output_Variables = NULL;
      Particle_Output_Variables = NULL;
      Particle_Filename_Template = NULL;
    }
    if (Particle_Number_PBCs)
      PBCs = (PBC_t *)calloc(Particle_Number_PBCs, sizeof(PBC_t));
    else
      PBCs = NULL;
  }

  return;
}

/*
 * noahs_ark() -- 2nd stage transport of init info from Proc 0 to world
 *
 * The bulk of the information scanned in is described to MPI in terms
 * of a monster derived data type.
 */

void noahs_ark(void) {
#ifdef PARALLEL
  int i;
  int j;
  int mode;
  MATRL_PROP_STRUCT *mp_ptr;
#endif

  DDD n = NULL; /* initialize to avoid picky complaints
                 * from compilers when PARALLEL is undefd */

#ifdef PARALLEL
  Noahs_Ark = ddd_alloc();
  n = Noahs_Ark;

  /*
   * Stuff read into rd_file_specs()
   */

  ddd_add_member(n, ExoFile, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, ExoFileOut, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, ExoFileOutMono, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, Init_GuessFile, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, Soln_OutFile, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, ExoAuxFile, MAX_FNL, MPI_CHAR);
  ddd_add_member(n, &ExoTimePlane, 1, MPI_INT);
  ddd_add_member(n, &Write_Intermediate_Solutions, 1, MPI_INT);
  ddd_add_member(n, &Write_Initial_Solution, 1, MPI_INT);

  if (GomaPetscOptionsStrLen > 0) {
    ddd_add_member(n, GomaPetscOptions, GomaPetscOptionsStrLen, MPI_CHAR);
  }
  /*
   * rd_genl_specs()
   */

  /*  ddd_add_member(n, &Num_Proc, 1, MPI_INT); ...dal */
  ddd_add_member(n, &Dim, 1, MPI_INT);
  ddd_add_member(n, &Iout, 1, MPI_INT);
  ddd_add_member(n, &Debug_Flag, 1, MPI_INT);
#ifdef MATRIX_DUMP
  ddd_add_member(n, &Number_Jac_Dump, 1, MPI_INT);
#endif
  ddd_add_member(n, &Guess_Flag, 1, MPI_INT);
  ddd_add_member(n, &Conformation_Flag, 1, MPI_INT);

  /*
   * The variable initialization structures are of fixed size, but only
   * communicate those that contain meaningful information read by proc 0
   * in mm_input.c
   */

  if (Num_Var_Init > 0 || Num_Var_Bound > 0) {
    for (i = 0; i < (Num_Var_Init + Num_Var_Bound); i++) {
      ddd_add_member(n, &Var_init[i].var, 1, MPI_INT);
      ddd_add_member(n, &Var_init[i].ktype, 1, MPI_INT);
      ddd_add_member(n, &Var_init[i].init_val, 1, MPI_DOUBLE);
      ddd_add_member(n, &Var_init[i].init_val_min, 1, MPI_DOUBLE);
      ddd_add_member(n, &Var_init[i].init_val_max, 1, MPI_DOUBLE);
    }
  }

  if (Num_Var_LS_Init > 0) {
    for (i = Num_Var_Init + Num_Var_Bound; i < Num_Var_Init + Num_Var_Bound + Num_Var_LS_Init;
         i++) {
      ddd_add_member(n, &Var_init[i].var, 1, MPI_INT);
      ddd_add_member(n, &Var_init[i].ktype, 1, MPI_INT);
      ddd_add_member(n, &Var_init[i].init_val_minus, 1, MPI_DOUBLE);
      ddd_add_member(n, &Var_init[i].init_val_plus, 1, MPI_DOUBLE);
    }
  }

  /*
   * Again, since the raven told us the actual number of external field
   * variables, we'll only communicate those.
   */

  if (efv->ev != T_NOTHING) {
    ddd_add_member(n, &efv->TALE, 1, MPI_INT);
    for (i = 0; i < efv->Num_external_field; i++) {
      ddd_add_member(n, efv->name[i], 20, MPI_CHAR);
      ddd_add_member(n, efv->file_nm[i], MAX_FNL, MPI_CHAR);
      ddd_add_member(n, &efv->i[i], 1, MPI_INT);
      ddd_add_member(n, efv->field_type[i], 15, MPI_CHAR);
    }
  }
  ddd_add_member(n, &efv->ev_etch_area, 1, MPI_INT);
  ddd_add_member(n, &efv->ev_etch_depth, 1, MPI_INT);

  ddd_add_member(n, &Anneal_Mesh, 1, MPI_INT);

  /*
   * rd_timeint_specs()
   */

  ddd_add_member(n, &tran->MaxTimeSteps, 1, MPI_INT);
  ddd_add_member(n, &tran->MaxSteadyStateSteps, 1, MPI_INT);
  ddd_add_member(n, &tran->Fill_Weight_Fcn, 1, MPI_INT);
  ddd_add_member(n, &tran->Fill_Equation, 1, MPI_INT);
  ddd_add_member(n, &tran->Delta_t0, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->Delta_t_min, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->Delta_t_max, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->TimeMax, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->theta, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->eps, 1, MPI_DOUBLE);

  /*
    for ( i=0; i<MAX_VARIABLE_TYPES; i++)
  */
  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    ddd_add_member(n, &tran->use_var_norm[i], 1, MPI_INT);
  }

  ddd_add_member(n, &tran->fix_freq, 1, MPI_INT);
  ddd_add_member(n, &tran->print_freq, 1, MPI_INT);
  ddd_add_member(n, &tran->print_delt, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->print_delt2_time, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->print_delt2, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->init_time, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->const_dt_after_failure, 1, MPI_INT);
  ddd_add_member(n, &tran->time_step_decelerator, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->resolved_delta_t_min, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->Courant_Limit, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->Restart_Time_Integ_After_Renorm, 1, MPI_INT);
  ddd_add_member(n, &tran->steady_state_tolerance, 1, MPI_DOUBLE);
  ddd_add_member(n, &tran->march_to_steady_state, 1, MPI_INT);
  ddd_add_member(n, &tran->ale_adapt, 1, MPI_INT);
  ddd_add_member(n, &tran->ale_adapt_freq, 1, MPI_INT);
  ddd_add_member(n, &tran->ale_adapt_iso_size, 1, MPI_DOUBLE);

  /*
   * Solver stuff
   */

  ddd_add_member(n, Matrix_Solver, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Format, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Scaling, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Preconditioner, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Residual_Norm_Type, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Output_Type, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Factorization_Reuse, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Maximum_Iterations, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Polynomial_Order, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Factor_Overlap, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Krylov_Subspace, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Orthogonalization, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Auxiliary_Vector, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Convergence_Tolerance, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Drop_Tolerance, MAX_CHAR_IN_INPUT, MPI_CHAR);

  /*
   * New for Aztec 2.0...
   */

  ddd_add_member(n, Matrix_Reorder, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Subdomain_Solver, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Graph_Fillin, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Overlap_Type, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Factorization_Save, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_ILUT_Fill_Factor, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_RILU_Relax_Factor, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_BILU_Threshold, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Relative_Threshold, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Matrix_Absolute_Threshold, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Amesos_Package, MAX_CHAR_IN_INPUT, MPI_CHAR);
  ddd_add_member(n, Stratimikos_File, MAX_CHAR_IN_INPUT * MAX_NUM_MATRICES, MPI_CHAR);

  ddd_add_member(n, &Linear_Solver, 1, MPI_INT);

  ddd_add_member(n, &Max_Newton_Steps, 1, MPI_INT);
  ddd_add_member(n, &Guess_Flag, 1, MPI_INT);
  ddd_add_member(n, &Conformation_Flag, 1, MPI_INT);
  ddd_add_member(n, &damp_factor1, 1, MPI_DOUBLE);
  ddd_add_member(n, &damp_factor2, 1, MPI_DOUBLE);
  ddd_add_member(n, &damp_factor3, 1, MPI_DOUBLE);
  ddd_add_member(n, &var_damp, MAX_VARIABLE_TYPES, MPI_DOUBLE);
  ddd_add_member(n, &custom_tol1, 1, MPI_DOUBLE);
  ddd_add_member(n, &custom_tol2, 1, MPI_DOUBLE);
  ddd_add_member(n, &custom_tol3, 1, MPI_DOUBLE);
  ddd_add_member(n, &Newt_Jacobian_Reformation_stride, 1, MPI_INT);
  ddd_add_member(n, &Time_Jacobian_Reformation_stride, 1, MPI_INT);
  ddd_add_member(n, &modified_newton, 1, MPI_INT);
  ddd_add_member(n, &convergence_rate_tolerance, 1, MPI_DOUBLE);
  ddd_add_member(n, &modified_newt_norm_tol, 1, MPI_DOUBLE);
  ddd_add_member(n, Epsilon, MAX_NUM_MATRICES * 3, MPI_DOUBLE);
  ddd_add_member(n, &Solver_Output_Format, 1, MPI_INT);
  ddd_add_member(n, &Output_Variable_Stats, 1, MPI_INT);

  /*
   * Eigensolver inputs.
   */

  ddd_add_member(n, &eigen->Eigen_Algorithm, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_NEV_WANT, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Maximum_Iterations, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Maximum_Outer_Iterations, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Filter, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Krylov_Subspace, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Recycle, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Record_Modes, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Matrix_Output, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Solve_Freq, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Write_Freq, 1, MPI_INT);
  ddd_add_member(n, &eigen->Eigen_Tolerance, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_IV_Wt, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Shifts[0], 5, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Cayley_Sigma, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Cayley_Mu, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_SI_Tol_Param, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Relative_Tol, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Linear_Tol, 1, MPI_DOUBLE);
  ddd_add_member(n, &eigen->Eigen_Output_File, 85, MPI_CHAR);

  /*
   * LSA 3D of 2D inputs.
   */
  ddd_add_member(n, &LSA_3D_of_2D_pass, 1, MPI_INT);
  ddd_add_member(n, &LSA_3D_of_2D_wave_number, 1, MPI_DOUBLE);
  ddd_add_member(n, &LSA_current_wave_number, 1, MPI_INT);
  for (i = 0; i < LSA_number_wave_numbers; i++) {
    ddd_add_member(n, &(LSA_wave_numbers[i]), 1, MPI_DOUBLE);
  }

  /*
   * Well, transport this over for now, then, later update each individual
   * processors conception of NZeroes.
   */

  ddd_add_member(n, &NZeros, 1, MPI_INT);

  ddd_add_member(n, &GNZeros, 1, MPI_INT);
  ddd_add_member(n, &PSPG, 1, MPI_INT);
  ddd_add_member(n, &PSPP, 1, MPI_INT);
  ddd_add_member(n, &PS_scaling, 1, MPI_DOUBLE);
  ddd_add_member(n, &Cont_GLS, 1, MPI_INT);
  ddd_add_member(n, &Filter_Species, 1, MPI_INT);
  ddd_add_member(n, &Include_Visc_Sens, 1, MPI_INT);
  ddd_add_member(n, &Visc_Sens_Copy, 1, MPI_INT);
  ddd_add_member(n, &Visc_Sens_Factor, 1, MPI_INT);
  ddd_add_member(n, &c_min, 1, MPI_DOUBLE);
  ddd_add_member(n, &c_max, 1, MPI_DOUBLE);

  /*
   * Particle data.
   */
  ddd_add_member(n, &Particle_Model, 1, MPI_INT);
  ddd_add_member(n, Particle_Model_Data, MAX_PARTICLE_MODEL_DATA_VALUES, MPI_DOUBLE);
  ddd_add_member(n, &Particle_Max_Time_Steps, 1, MPI_INT);
  ddd_add_member(n, &Particle_Number, 1, MPI_INT);
  ddd_add_member(n, Particle_Restart_Filename, MAX_PARTICLE_FILENAME_LENGTH, MPI_CHAR);
  ddd_add_member(n, &Particle_Output_Stride, 1, MPI_INT);
  ddd_add_member(n, &Particle_Output_Time_Step, 1, MPI_DOUBLE);
  ddd_add_member(n, &Particle_Output_Format, 1, MPI_INT);
  ddd_add_member(n, &Particle_Density, 1, MPI_DOUBLE);
  ddd_add_member(n, &Particle_Radius, 1, MPI_DOUBLE);
  ddd_add_member(n, &Particle_Ratio, 1, MPI_DOUBLE);
  ddd_add_member(n, &Particle_Creation_Domain, 1, MPI_INT);
  ddd_add_member(n, &Particle_Move_Domain, 1, MPI_INT);
  ddd_add_member(n, Particle_Creation_Domain_Reals, MAX_DOMAIN_REAL_VALUES, MPI_DOUBLE);
  ddd_add_member(n, Particle_Move_Domain_Reals, MAX_DOMAIN_REAL_VALUES, MPI_DOUBLE);
  ddd_add_member(n, Particle_Creation_Domain_Filename, MAX_PARTICLE_FILENAME_LENGTH, MPI_CHAR);
  ddd_add_member(n, Particle_Creation_Domain_Name, MAX_PARTICLE_STRING_LENGTH, MPI_CHAR);
  ddd_add_member(n, Particle_Move_Domain_Filename, MAX_PARTICLE_FILENAME_LENGTH, MPI_CHAR);
  ddd_add_member(n, Particle_Move_Domain_Name, MAX_PARTICLE_STRING_LENGTH, MPI_CHAR);
  ddd_add_member(n, &Particle_Full_Output_Stride, 1, MPI_INT);
  ddd_add_member(n, Particle_Full_Output_Filename, MAX_PARTICLE_FILENAME_LENGTH, MPI_CHAR);
  if (Particle_Number_Sample_Types) {
    ddd_add_member(n, Particle_Number_Samples, Particle_Number_Sample_Types, MPI_INT);
    ddd_add_member(n, Particle_Number_Output_Variables, Particle_Number_Sample_Types, MPI_INT);
    for (i = 0; i < Particle_Number_Sample_Types; i++)
      ddd_add_member(n, Particle_Filename_Template[i], MAX_PARTICLE_FILENAME_LENGTH, MPI_CHAR);
  }
  for (i = 0; i < Particle_Number_PBCs; i++) {
    ddd_add_member(n, &PBCs[i].SS_id, 1, MPI_INT);
    ddd_add_member(n, &PBCs[i].type, 1, MPI_INT);
    ddd_add_member(n, PBCs[i].real_data, MAX_PBC_REAL_VALUES, MPI_DOUBLE);
    ddd_add_member(n, PBCs[i].int_data, MAX_PBC_INT_VALUES, MPI_INT);
    ddd_add_member(n, PBCs[i].string_data, MAX_PBC_STRING_DATA_LENGTH, MPI_CHAR);
  }

  /*
   * rd_geometry_specs line ... No, really, it was smarter to do it this way.
   */

  /*
   * Boundary Conditions.
   */

  /*
   * Set up skeleton of the BC_Types...
   */

  for (i = 0; i < Num_BC; i++) {
    ddd_add_member(n, &BC_Types[i].BC_Name, 1, MPI_INT);
    ddd_add_member(n, BC_Types[i].Set_Type, 3, MPI_CHAR);
    ddd_add_member(n, &BC_Types[i].BC_ID, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_ID2, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_ID3, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_matrl_index_1, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_matrl_index_2, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_matrl_index_3, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_matrl_index_4, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].BC_EBID_Apply, 1, MPI_INT);
    ddd_add_member(n, BC_Types[i].BC_Data_Int, MAX_BC_INT_DATA, MPI_INT);
    ddd_add_member(n, BC_Types[i].BC_Data_Float, MAX_BC_FLOAT_DATA, MPI_DOUBLE);
    ddd_add_member(n, &BC_Types[i].len_u_BC, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].max_DFlt, 1, MPI_INT);

    /* double u_BC delivered by Dove...*/

    ddd_add_member(n, &BC_Types[i].BC_Desc_index, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].index_dad, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].equation, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].matrix, 1, MPI_INT);
    ddd_add_member(n, &BC_Types[i].species_eq, 1, MPI_INT);

    ddd_add_member(n, &BC_Types[i].BC_relax, 1, MPI_DOUBLE);
    ddd_add_member(n, &BC_Types[i].table_index, 1, MPI_INT);
  }

  /*
   * This will transport all of the new dynamically allocated BC Descriptions
   * to other processors *except* for the char *name1, *name2.
   *
   * Currently, those names are pointers to other character arrays in the
   * statically declared bundle of BC Descriptions. That's great. We'll just
   * use the BC_Desc_index to find out which one so that these pointers
   * can be reconstructed on the other processors. This is done
   * after the ark is sent in the "ark_landing" below.
   */

  for (i = 0; i < num_new_BC_Desc; i++) {
    ddd_add_member(n, &(new_BC_Desc[i]->method), 1, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->BC_Name), 1, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->equation), 1, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->vector), 1, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->rotate), 1, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->sens[0]), MAX_VARIABLE_TYPES, MPI_INT);
    ddd_add_member(n, &(new_BC_Desc[i]->i_apply), 1, MPI_INT);
  }

  /*
   * Now we have to include the Table data.  We dont carry along the
   * names for ordinates and abscissas to other processors.
   * The table data itself will be allocated in "ark_landing"
   */

  for (i = 0; i < num_AC_Tables; i++) {
    ddd_add_member(n, AC_Tables[i]->t_name[0], 132, MPI_CHAR);
    ddd_add_member(n, AC_Tables[i]->t_name[1], 132, MPI_CHAR);
    ddd_add_member(n, AC_Tables[i]->t_index, 2, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->columns), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->interp_method), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->tablelength), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->f_index), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->species_eq), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->ngrid), 1, MPI_INT);
    ddd_add_member(n, &(AC_Tables[i]->yscale), 1, MPI_DOUBLE);
    ddd_add_member(n, &(AC_Tables[i]->Emin), 1, MPI_DOUBLE);
  }

  for (i = 0; i < num_BC_Tables; i++) {
    ddd_add_member(n, BC_Tables[i]->t_name[0], 132, MPI_CHAR);
    ddd_add_member(n, BC_Tables[i]->t_name[1], 132, MPI_CHAR);
    ddd_add_member(n, BC_Tables[i]->t_index, 2, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->columns), 1, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->interp_method), 1, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->tablelength), 1, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->f_index), 1, MPI_INT);
    /*    BC_Tables[i]->slope is computed locally  */
    ddd_add_member(n, &(BC_Tables[i]->species_eq), 1, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->ngrid), 1, MPI_INT);
    ddd_add_member(n, &(BC_Tables[i]->yscale), 1, MPI_DOUBLE);
    ddd_add_member(n, &(BC_Tables[i]->Emin), 1, MPI_DOUBLE);
  }

  for (i = 0; i < num_MP_Tables; i++) {
    ddd_add_member(n, MP_Tables[i]->t_name[0], 132, MPI_CHAR);
    ddd_add_member(n, MP_Tables[i]->t_name[1], 132, MPI_CHAR);
    ddd_add_member(n, MP_Tables[i]->t_index, 2, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->columns), 1, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->interp_method), 1, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->tablelength), 1, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->f_index), 1, MPI_INT);
    /*    MP_Tables[i]->slope is computed locally    */
    ddd_add_member(n, &(MP_Tables[i]->species_eq), 1, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->ngrid), 1, MPI_INT);
    ddd_add_member(n, &(MP_Tables[i]->yscale), 1, MPI_DOUBLE);
    ddd_add_member(n, &(MP_Tables[i]->Emin), 1, MPI_DOUBLE);
  }

  for (i = 0; i < num_ext_Tables; i++) {
    ddd_add_member(n, ext_Tables[i]->t_name[0], 132, MPI_CHAR);
    ddd_add_member(n, ext_Tables[i]->t_name[1], 132, MPI_CHAR);
    ddd_add_member(n, ext_Tables[i]->t_name[2], 132, MPI_CHAR);
    ddd_add_member(n, ext_Tables[i]->t_index, 3, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->columns), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->interp_method), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->tablelength), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->f_index), 1, MPI_INT);
    /*    ext_Tables[i]->slope is computed locally  */
    ddd_add_member(n, &(ext_Tables[i]->species_eq), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->ngrid), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->ngrid2), 1, MPI_INT);
    ddd_add_member(n, &(ext_Tables[i]->yscale), 1, MPI_DOUBLE);
    ddd_add_member(n, &(ext_Tables[i]->Emin), 1, MPI_DOUBLE);
  }

  ddd_add_member(n, &PRESSURE_DATUM, 1, MPI_INT);
  ddd_add_member(n, &pressure_datum_element, 1, MPI_INT);
  ddd_add_member(n, &pressure_datum_value, 1, MPI_DOUBLE);

  /*
   * Rotation Conditions. Initially assume that Rotation Specifications
   * contain BC_Desc[] references to statically declared data. Thus, a
   * simple integer index suffices to re-create the information on other
   * processors.
   *
   * If you need to reference new, different dynamically allocated versions
   * of BC_Desc[], then see the BC_Types[] conniptions above for a method
   * of implementing this.
   */

  for (i = 0; i < Num_ROT; i++) {
    ddd_add_member(n, &ROT_Types[i].eq_type, 1, MPI_INT);
    ddd_add_member(n, &ROT_Types[i].type, 1, MPI_INT);
    ddd_add_member(n, ROT_Types[i].ss_id, CDIM, MPI_INT);
    ddd_add_member(n, ROT_Types[i].BC_Type, CDIM, MPI_INT);
    ddd_add_member(n, ROT_Types[i].BC_SS, CDIM, MPI_INT);
    ddd_add_member(n, ROT_Types[i].BC_id, CDIM, MPI_INT);

    /*
     * Reconstruct the ptrs to BC_desc's via these integer indeces...
     */

    ddd_add_member(n, ROT_Types[i].BC_desc_index, CDIM, MPI_INT);

    ddd_add_member(n, &ROT_Types[i].method, 1, MPI_INT);
    ddd_add_member(n, &ROT_Types[i].ROTATE, 1, MPI_INT);
    ddd_add_member(n, ROT_Types[i].seed, CDIM, MPI_DOUBLE);
    ddd_add_member(n, &ROT_Types[i].node, 1, MPI_INT);

    /*
     * The ptr to a list of elements on the side set is empty for now.
     * It gets filled later in bc/rotate.c; when it does it will mean
     * different things to different processors...
     */
    /*
          ddd_add_member(n, &ROT_Types[i].num_elem, 1, MPI_INT);
          ddd_add_member(n, &ROT_Types[i].ss_ptr, 1, MPI_INT);
    */
  }

  /*
   * Materials and Equation Specifications...
   */
  /*
   * Broadcast the Uniform_Problem_Description Information
   */
  ddd_add_member(n, &upd->Total_Num_EQ, MAX_NUM_MATRICES, MPI_INT);
  ddd_add_member(n, &upd->Total_Num_Var, MAX_NUM_MATRICES, MPI_INT);
  ddd_add_member(n, &upd->CoordinateSystem, 1, MPI_INT);
  ddd_add_member(n, upd->vp, MAX_NUM_MATRICES * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_INT);
  ddd_add_member(n, upd->ep, MAX_NUM_MATRICES * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_INT);
  ddd_add_member(n, &upd->Max_Num_Species, 1, MPI_INT);
  ddd_add_member(n, &upd->Max_Num_Species_Eqn, 1, MPI_INT);
  ddd_add_member(n, &upd->Tot_Num_VolSpecies, 1, MPI_INT);
  ddd_add_member(n, &upd->Num_Mat, 1, MPI_INT);
  ddd_add_member(n, &upd->Species_Var_Type, 1, MPI_INT);
  ddd_add_member(n, &upd->Pressure_Datum, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->Max_Num_Porous_Eqn, 1, MPI_INT);
  ddd_add_member(n, &upd->XFEM, 1, MPI_INT);
  ddd_add_member(n, &upd->Process_Temperature, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->Acoustic_Frequency, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->EM_Frequency, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->Free_Space_Permittivity, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->Free_Space_Permeability, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->Light_Cosmu, 1, MPI_DOUBLE);
  ddd_add_member(n, &upd->SegregatedSolve, 1, MPI_INT);
  ddd_add_member(n, &upd->SegregatedSubcycles, 1, MPI_INT);
  ddd_add_member(n, &upd->PSPG_advection_correction, 1, MPI_INT);
  ddd_add_member(n, &upd->petsc_solve_post_proc, 1, MPI_INT);
  ddd_add_member(n, &upd->devss_traceless_gradient, 1, MPI_INT);

  ddd_add_member(n, pg->time_step_control_disabled, MAX_NUM_MATRICES, MPI_INT);
  ddd_add_member(n, pg->matrix_subcycle_count, MAX_NUM_MATRICES, MPI_INT);
  ddd_add_member(n, upd->matrix_index, MAX_EQNS, MPI_INT);

  // turbulence information
  ddd_add_member(n, &(upd->turbulent_info->use_internal_wall_distance), 1, MPI_INT);
  ddd_add_member(n, &(upd->turbulent_info->num_node_sets), 1, MPI_INT);
  ddd_add_member(n, &(upd->turbulent_info->num_side_sets), 1, MPI_INT);

  for (i = 0; i < upd->Num_Mat; i++) {
    int imtrx;
    ddd_add_member(n, &pd_glob[i]->Num_EQ, MAX_NUM_MATRICES, MPI_INT);
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      ddd_add_member(n, pd_glob[i]->e[imtrx], MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->v[imtrx], MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->w[imtrx], MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->i[imtrx], MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->m[imtrx], MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->gv, MAX_EQNS, MPI_INT);
      ddd_add_member(n, pd_glob[i]->mi, MAX_EQNS, MPI_INT);

      ddd_add_member(n, &pd_glob[i]->etm[imtrx][0][0], MAX_EQNS * MAX_TERM_TYPES, MPI_DOUBLE);
    }
    /*  pd_glob[i]->CoordinateSystem  is already assigned */
    ddd_add_member(n, &pd_glob[i]->MeshMotion, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->MeshInertia, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->RealSolidFluxModel, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->MassFluxModel, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->MomentumFluxModel, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->PorousFluxModel, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Num_Dim, 1, MPI_INT);
    /* pd_glob[i]->TimeIntegration  is already assigned */
    ddd_add_member(n, &pd_glob[i]->Continuation, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->AugmentingConditions, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->IntegrationMap, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->ShapeVar, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->ProjectionVar, 1, MPI_INT);
    /* pd_glob[i]->Num_Mat          is already assigned */
    ddd_add_member(n, pd_glob[i]->MaterialName, MAX_MATLNAME, MPI_CHAR);
    ddd_add_member(n, &pd_glob[i]->Num_Species, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Num_Species_Eqn, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Species_Var_Type, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Num_Rxn, 1, MPI_INT);

    /* material flags associated with augmented conditions */

    ddd_add_member(n, &pd_glob[i]->VolumeIntegral, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->LSVelocityIntegral, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Num_Porous_Eqn, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Num_Porous_Shell_Eqn, 1, MPI_INT);
    ddd_add_member(n, &pd_glob[i]->Do_Surf_Geometry, 1, MPI_INT);
  }

  /*
   * material block Initialization info
   */

  ddd_add_member(n, Num_Var_Init_Mat, MAX_NUMBER_MATLS, MPI_INT);
  for (i = 0; i < MAX_NUMBER_MATLS; i++) {
    for (j = 0; j < MAX_VARIABLE_TYPES + MAX_CONC; j++) {
      ddd_add_member(n, &Var_init_mat[i][j].var, 1, MPI_INT);
      ddd_add_member(n, &Var_init_mat[i][j].ktype, 1, MPI_INT);
      ddd_add_member(n, &Var_init_mat[i][j].init_val, 1, MPI_DOUBLE);
      ddd_add_member(n, &Var_init_mat[i][j].slave_block, 1, MPI_INT);
      ddd_add_member(n, &Var_init_mat[i][j].len_u_pars, 1, MPI_INT);
    }
  }

  /*
   * Include the Element quality metric data.
   * This structure was previously allocated in "raven_landing"
   */

  if (nEQM > 0) {
    ddd_add_member(n, &eqm->do_jac, 1, MPI_INT);
    ddd_add_member(n, &eqm->do_vol, 1, MPI_INT);
    ddd_add_member(n, &eqm->do_ang, 1, MPI_INT);
    ddd_add_member(n, &eqm->do_tri, 1, MPI_INT);
    ddd_add_member(n, &eqm->wt_jac, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->wt_vol, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->wt_ang, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->wt_tri, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_jac, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_vol, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_ang, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_tri, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_avg, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_low, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->eq_tol, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->vol_sum, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->vol_low, 1, MPI_DOUBLE);
    ddd_add_member(n, &eqm->vol_count, 1, MPI_INT);
    ddd_add_member(n, &eqm->tol_type, 1, MPI_INT);
  }

  /*
   * Include the Augmenting condition data.
   * This structure was previously allocated in "raven_landing"
   */

  for (i = 0; i < nAC; i++) {
    ddd_add_member(n, &augc[i].nAC, 1, MPI_INT);
    ddd_add_member(n, &augc[i].iread, 1, MPI_INT);
    ddd_add_member(n, &augc[i].theta, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].eps, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].Type, 1, MPI_INT);
    ddd_add_member(n, &augc[i].BCID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].DFID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].DHID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].VOLID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].COMPID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].SSID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].MTID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].MPID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].MDID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].MFID, 1, MPI_INT);
    ddd_add_member(n, &augc[i].len_AC, 1, MPI_INT);
    /*   DataFlt allocated in ark_landing */
    ddd_add_member(n, &augc[i].tmp1, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].tmp2, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].tmp3, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].evol, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].CONSTV, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].LewisNum, 1, MPI_DOUBLE);
    /*   derivative block d_evol_dx is later allocated in rf_solve */
    ddd_add_member(n, &augc[i].LSPHASE, 1, MPI_INT);
    ddd_add_member(n, &augc[i].DIR, 1, MPI_INT);
    ddd_add_member(n, &augc[i].fluid_eb, 1, MPI_INT);
    ddd_add_member(n, &augc[i].solid_eb, 1, MPI_INT);
    ddd_add_member(n, &augc[i].lm_eb, 1, MPI_INT);
    ddd_add_member(n, &augc[i].lm_elem, 1, MPI_INT);
    ddd_add_member(n, &augc[i].lm_side, 1, MPI_INT);
    ddd_add_member(n, &augc[i].lm_dim, 1, MPI_INT);
    ddd_add_member(n, &augc[i].lm_value, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].lm_resid, 1, MPI_DOUBLE);
    ddd_add_member(n, &augc[i].Params_File, 128, MPI_CHAR);
    ddd_add_member(n, &augc[i].AP_param, 64, MPI_CHAR);
    ddd_add_member(n, &augc[i].Aprepro_lib_string_len, 1, MPI_INT);
  }

  /*
   * Include the Continuation information
   */

  if (Continuation != ALC_NONE) {
    ddd_add_member(n, &cont->MaxPathSteps, 1, MPI_INT);
    ddd_add_member(n, &cont->PathIntr, 1, MPI_INT);
    ddd_add_member(n, &cont->Delta_s0, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->Delta_s_min, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->Delta_s_max, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->PathMax, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->alpha, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->beta, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->gamma, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->delta, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->theta, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->eps, 1, MPI_DOUBLE);
    ddd_add_member(n, cont->use_var_norm, MAX_VARIABLE_TYPES, MPI_INT);
    ddd_add_member(n, &cont->print_freq, 1, MPI_INT);
    ddd_add_member(n, &cont->fix_freq, 1, MPI_INT);
    ddd_add_member(n, &cont->print_delt, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->print_delt2_path, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->print_delt2, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->radius, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->BegParameterValue, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->EndParameterValue, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->InitDir, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->upType, 1, MPI_INT);
    ddd_add_member(n, &cont->upBCID, 1, MPI_INT);
    ddd_add_member(n, &cont->upDFID, 1, MPI_INT);
    ddd_add_member(n, &cont->upDHID, 1, MPI_INT);
    ddd_add_member(n, &cont->upMTID, 1, MPI_INT);
    ddd_add_member(n, &cont->upMPID, 1, MPI_INT);
    ddd_add_member(n, &cont->upMDID, 1, MPI_INT);
    ddd_add_member(n, &cont->upMFID, 1, MPI_INT);
    ddd_add_member(n, &cont->sensvec_id, 1, MPI_INT);
    ddd_add_member(n, &cont->tmp1, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->tmp2, 1, MPI_DOUBLE);
    ddd_add_member(n, &cont->tmp3, 1, MPI_DOUBLE);
  }

  /*
   * Include LOCA inputs
   */

  /* This one may apply even without continuation! */
  ddd_add_member(n, &loca_in->Cont_Alg, 1, MPI_INT);

  if (Continuation == LOCA) {
    ddd_add_member(n, &loca_in->Cont_Order, 1, MPI_INT);
    ddd_add_member(n, &loca_in->StepAggr, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->debug, 1, MPI_INT);
    ddd_add_member(n, &loca_in->DpDs2, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->DpDsHi, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->Texp, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->MaxTS, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->perturb, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->TPupType, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPupBCID, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPupDFID, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPupMTID, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPupMPID, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPupMDID, 1, MPI_INT);
    ddd_add_member(n, &loca_in->TPGuess, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->TPFinal, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->NVRestart, 1, MPI_INT);
    ddd_add_member(n, &loca_in->NV_exoII_infile, MAX_FNL, MPI_CHAR);
    ddd_add_member(n, &loca_in->NV_imag_infile, MAX_FNL, MPI_CHAR);
    ddd_add_member(n, &loca_in->omega, 1, MPI_DOUBLE);
    ddd_add_member(n, &loca_in->NVSave, 1, MPI_INT);
    ddd_add_member(n, &loca_in->NV_exoII_outfile, MAX_FNL, MPI_CHAR);
    ddd_add_member(n, &loca_in->NV_imag_outfile, MAX_FNL, MPI_CHAR);
    ddd_add_member(n, &loca_in->NV_time_index, 1, MPI_INT);

    /*
     * Multiple continuation conditions
     */

    for (i = 0; i < nCC; i++) {
      ddd_add_member(n, &cpcc[i].nCC, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].ratio, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].old_value, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].value, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].Type, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].fn_flag, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].BCID, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].DFID, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].MTID, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].MPID, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].MDID, 1, MPI_INT);
      ddd_add_member(n, &cpcc[i].Beg_CC_Value, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].End_CC_Value, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].coeff_0, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].coeff_1, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].coeff_2, 1, MPI_DOUBLE);
      ddd_add_member(n, &cpcc[i].sensvec_id, 1, MPI_INT);
    }

    for (i = 0; i < nTC; i++) {
      ddd_add_member(n, &tpcc[i].nCC, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].ratio, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].old_value, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].value, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].Type, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].fn_flag, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].BCID, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].DFID, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].MTID, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].MPID, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].MDID, 1, MPI_INT);
      ddd_add_member(n, &tpcc[i].Beg_CC_Value, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].End_CC_Value, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].coeff_0, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].coeff_1, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].coeff_2, 1, MPI_DOUBLE);
      ddd_add_member(n, &tpcc[i].sensvec_id, 1, MPI_INT);
    }

  } /* End if (Continuation == LOCA) */

  /*
   * Include the Hunting condition data.
   * This structure was previously allocated in "raven_landing"
   */

  for (i = 0; i < nHC; i++) {
    ddd_add_member(n, &hunt[i].nHC, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].theta, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].eps, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].Type, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].ramp, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].BCID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].DFID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].DHID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].MTID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].MPID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].MDID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].MFID, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].BegParameterValue, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].EndParameterValue, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].Delta_s0, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].Delta_s_min, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].Delta_s_max, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].sensvec_id, 1, MPI_INT);
    ddd_add_member(n, &hunt[i].tmp1, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].tmp2, 1, MPI_DOUBLE);
    ddd_add_member(n, &hunt[i].tmp3, 1, MPI_DOUBLE);
  }

  /*
   * Just transport the interpolations that this problem will actually use...
   */

  for (i = 0; i < Num_Interpolations; i++) {
    ddd_add_member(n, Unique_Interpolations, Num_Interpolations, MPI_INT);
  }

  if (ls != NULL) {
    ddd_add_member(n, &ls->var, 1, MPI_INT);
    ddd_add_member(n, &ls->Use_Level_Set, 1, MPI_INT);
    ddd_add_member(n, &ls->Evolution, 1, MPI_INT);
    ddd_add_member(n, &ls->Contact_Inflection, 1, MPI_INT);
    ddd_add_member(n, &ls->Isosurface_Subsurf_Type, 1, MPI_INT);
    ddd_add_member(n, &ls->Init_Method, 1, MPI_INT);
    ddd_add_member(n, &ls->Num_Var_Init, 1, MPI_INT);
    ddd_add_member(n, &ls->Length_Scale, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->adapt, 1, MPI_INT);
    ddd_add_member(n, &ls->adapt_freq, 1, MPI_INT);
    ddd_add_member(n, &ls->adapt_inner_size, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->adapt_outer_size, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->adapt_width, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Control_Width, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Renorm_Tolerance, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Renorm_Method, 1, MPI_INT);
    ddd_add_member(n, &ls->Search_Option, 1, MPI_INT);
    ddd_add_member(n, &ls->Grid_Search_Depth, 1, MPI_INT);
    ddd_add_member(n, &ls->Integration_Depth, 1, MPI_INT);
    ddd_add_member(n, &ls->Interface_Output, 1, MPI_INT);
    ddd_add_member(n, &ls->Renorm_Freq, 1, MPI_INT);
    ddd_add_member(n, &ls->Renorm_Countdown, 1, MPI_INT);
    ddd_add_member(n, &ls->Force_Initial_Renorm, 1, MPI_INT);
    ddd_add_member(n, &ls->Initial_LS_Displacement, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Mass_Value, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Mass_Sign, 1, MPI_INT);
    ddd_add_member(n, &ls->Contact_Tolerance, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Fluid_Solid, 1, MPI_INT);
    ddd_add_member(n, &ls->Fluid_Sign, 1, MPI_INT);
    ddd_add_member(n, &ls->Solid_Sign, 1, MPI_INT);
    ddd_add_member(n, &ls->Elem_Sign, 1, MPI_INT);
    ddd_add_member(n, &ls->SubElemIntegration, 1, MPI_INT);
    ddd_add_member(n, &ls->AdaptIntegration, 1, MPI_INT);
    ddd_add_member(n, &ls->Adaptive_Order, 1, MPI_INT);
    ddd_add_member(n, &ls->Ghost_Integ, 1, MPI_INT);
    ddd_add_member(n, &ls->Ghost_Integ_Active, 1, MPI_INT);
    ddd_add_member(n, &ls->CrossMeshQuadPoints, 1, MPI_INT);
    ddd_add_member(n, &ls->Extension_Velocity, 1, MPI_INT);
    ddd_add_member(n, &ls->CalcSurfDependencies, 1, MPI_INT);
    ddd_add_member(n, &ls->Ignore_F_deps, 1, MPI_INT);
    ddd_add_member(n, &ls->Periodic_Planes, 1, MPI_INT);
    ddd_add_member(n, &ls->Periodic_Plane_Loc, 6, MPI_DOUBLE);
    ddd_add_member(n, &ls->PSPP_filter, 1, MPI_INT);
    ddd_add_member(n, &ls->Sat_Hyst_Renorm_Lockout, 1, MPI_INT);
    ddd_add_member(n, &ls->SubcyclesAfterRenorm, 1, MPI_INT);
    ddd_add_member(n, &ls->ghost_stress, 1, MPI_INT);
    ddd_add_member(n, &ls->Toure_Penalty, 1, MPI_INT);
    ddd_add_member(n, &ls->YZbeta, 1, MPI_INT);
    ddd_add_member(n, &ls->YZbeta_scale, 1, MPI_DOUBLE);
    ddd_add_member(n, &ls->Huygens_Freeze_Nodes, 1, MPI_INT);
    ddd_add_member(n, &ls->Semi_Implicit_Integration, 1, MPI_INT);
  }

  if (pfd != NULL) {
    ddd_add_member(n, &pfd->num_phase_funcs, 1, MPI_INT);
    ddd_add_member(n, &pfd->Use_Phase_Field, 1, MPI_INT);
    for (i = 0; i < pfd->num_phase_funcs; i++) {
      ddd_add_member(n, &pfd->ls[i]->var, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Use_Level_Set, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Evolution, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Contact_Inflection, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Isosurface_Subsurf_Type, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Init_Method, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Num_Var_Init, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Length_Scale, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Control_Width, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Renorm_Tolerance, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Renorm_Method, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Search_Option, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Grid_Search_Depth, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Integration_Depth, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Interface_Output, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Renorm_Freq, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Renorm_Countdown, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Force_Initial_Renorm, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Initial_LS_Displacement, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Mass_Value, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Mass_Sign, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Contact_Tolerance, 1, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->Fluid_Solid, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Fluid_Sign, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Solid_Sign, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Elem_Sign, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->SubElemIntegration, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->AdaptIntegration, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Adaptive_Order, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Ghost_Integ, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Ghost_Integ_Active, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->CrossMeshQuadPoints, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Extension_Velocity, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->CalcSurfDependencies, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Ignore_F_deps, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Periodic_Planes, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Periodic_Plane_Loc, 6, MPI_DOUBLE);
      ddd_add_member(n, &pfd->ls[i]->PSPP_filter, 1, MPI_INT);
      ddd_add_member(n, &pfd->ls[i]->Sat_Hyst_Renorm_Lockout, 1, MPI_INT);
    }
  }

  /*
   * Material properties...
   */

  for (i = 0; i < upd->Num_Mat; i++) {
    mp_ptr = mp_glob[i];
    /*
     * Scalar property values, Booleans, flags, etc.
     */
    ddd_add_member(n, &mp_ptr->MatID, 1, MPI_INT);
    ddd_add_member(n, mp_ptr->Material_Name, MAX_MATLNAME, MPI_CHAR);
    ddd_add_member(n, &mp_ptr->Num_Matrl_Elem_Blk, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->DefaultDatabase, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->Num_Species, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->Num_Species_Eqn, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->Dropped_Last_Species_Eqn, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->NonDiluteFormulation, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->Num_Porous_Eqn, 1, MPI_INT);
    ddd_add_member(n, &mp_ptr->Num_Porous_Shell_Eqn, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->FlowingLiquid_viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Inertia_coefficient, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Spwt_func, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->SpSSPG_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SpSSPG_func, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->SpYZbeta_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SpYZbeta_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->SpYZbeta_value, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->Porous_wt_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Volume_Expansion, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->current_source, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->Momentwt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Momentwt_func, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->MomentSSPG_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MomentSSPG_func, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->MomentShock_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MomentShock_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->MomentShock_Ref, MAX_MOMENTS, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->MomentSecondLevelSetDiffusivityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MomentSecondLevelSetDiffusivity, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->MomentLevelSetDiffusionOnly, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->moment_source, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->density, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->electrical_conductivity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->permittivity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->permittivity_imag, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->magnetic_permeability, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->incident_wave, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->VoltageFormulation, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->heat_capacity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->heat_source, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->heat_transfer_coefficient, 1, MPI_DOUBLE);
    /*
     * The elc_rs_glob[i] evidently transport these Lame coefficients now.
     *
     * ddd_add_member(n, &mp_glob[i]->lame_lambda, 1, MPI_DOUBLE);
     * ddd_add_member(n, &mp_glob[i]->lame_mu, 1, MPI_DOUBLE);
     */
    ddd_add_member(n, &mp_glob[i]->mass_source, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->melting_point_liquidus, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->melting_point_solidus, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->porosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->porous_compressibility, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->initial_porosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->matrix_density, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->specific_heat, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->permeability, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->permeability_imag, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousLiqCompress, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousLiqRefPress, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->rel_gas_perm, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->rel_liq_perm, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->saturation, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->cap_pres, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->surface_tension, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->SurfaceDiffusionCoeffProjectionEqn, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->thermal_conductivity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Ewt_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Rst_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Rst_diffusion, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Rst_func_supg, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Mwt_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->SAwt_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->dilationalViscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->dilationalViscosityRatio, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->dilationalViscosityMultiplier, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->volumeFractionGas, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->reaction_rate, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->solution_temperature, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->electrolyte_temperature, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->electrolyte_conductivity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->electrode_potential, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->thermodynamic_potential, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->interfacial_area, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedPorosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedRadius, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedHeight, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedP0, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellPatm, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellPref, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellCrossKappa, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellInitPorePres, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellDiffusivity, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellRT, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellHenry, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->heightU, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->heightL, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->veloU[0], DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->veloL[0], DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->dcaU, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->dcaL, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->FilmEvap, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->DisjPress, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->DiffCoeff, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->lubsource, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->lubmomsource[0], DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->shell_user_par, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->acoustic_impedance, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->wave_number, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->acoustic_absorption, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->acoustic_ksquared_sign, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->refractive_index, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->light_absorption, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->extinction_index, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->moment_coalescence_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->moment_coalescence_scale, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->moment_growth_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->moment_growth_scale, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->moment_growth_reference_pressure, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->CapStress, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ConductivityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Ewt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Energy_Div_Term, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Rst_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Mwt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SAwt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->CurrentSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MomentSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->HeightUFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->HeightLFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->VeloUFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->VeloLFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DcaUFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DcaLFunctionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->FSIModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->TurbulentLubricationModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LubIntegrationModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedPorosityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedRadiusModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedHeightModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellClosedP0Model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellPatmModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellPrefModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellCrossKappaModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellInitPorePresModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellDiffusivityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellRTModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousShellHenryModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->FilmEvapModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DisjPressModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DiffCoeffModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LubSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LubMomSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DensityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->table_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Elec_ConductivityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->FlowingLiquidViscosityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->HeatCapacityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->HeatSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->InertiaCoefficientModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LiquidusModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MassSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousMediaType, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MeshSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->RealSolidSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MomentumSourceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousCompressibilityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PermeabilityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousLiquidCompressModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousLiqRefPressModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorosityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousMatrixDensityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousSpecificHeatModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->RelGasPermModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->RelLiqPermModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SaturationModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->CapPresModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Porous_wt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Porous_Mass_Lump, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SolidusModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Spwt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ReactionRateModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SolutionTemperatureModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SurfaceTensionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ViscosityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->DilationalViscosityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->VolumeExpansionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->i_ys, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ThermodynamicPotentialModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->InterfacialAreaModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->wave_numberModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Acoustic_ImpedanceModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Acoustic_AbsorptionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Refractive_IndexModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Extinction_IndexModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Light_AbsorptionModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Shell_User_ParModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PermittivityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->MagneticPermeabilityModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->IncidentWaveModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PBE_BA_Type, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SBM_Length_enabled, 1, MPI_INT);

    /* External field indeces PRS 10-1-2013 (shutdown times) */

    ddd_add_member(n, &mp_glob[i]->porosity_external_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->perm_external_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Xperm_external_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->SAT_external_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->cap_pres_external_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->por_shell_closed_porosity_ext_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->por_shell_closed_height_ext_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->por_shell_closed_radius_ext_field_index, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->heightU_ext_field_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->heightL_ext_field_index, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->thermal_cond_external_field, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->elec_cond_external_field, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->species_source_external_field_index, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->dist_wall_ext_field_index, 1, MPI_INT);

    /*
     * Material properties that are fixed length vectors of doubles.
     * Lengths are: total variable types (regular + concentrations)
     */

    ddd_add_member(n, mp_glob[i]->d2_viscosity, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_FlowingLiquid_viscosity, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_Inertia_coefficient, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_Volume_Expansion, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_current_source, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_moment_source, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_density, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_electrical_conductivity, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_heat_capacity, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_heat_source, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_heat_transfer_coefficient, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_mass_source, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_melting_point_liquidus, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_melting_point_solidus, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_porosity, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_permeability, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_PorousLiquidCompres, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_PorousLiqRefPress, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_rel_gas_perm, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_rel_liq_perm, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_saturation, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_cap_pres, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_species_source, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_surface_tension, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_thermal_conductivity, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_viscosity, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_dilationalViscosityRatio, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_volumeFractionGas, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_reaction_rate, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_solution_temperature, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_thermodynamic_potential, MAX_VARIABLE_TYPES + MAX_CONC,
                   MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_shell_user_par, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_acoustic_impedance, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_wave_number, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_acoustic_absorption, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_refractive_index, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_light_absorption, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_extinction_index, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->d_permittivity, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    /*
     * Material properties that are fixed length matrices of doubles.
     * Lengths are: total variable types (regular + concentrations)**2
     */

    ddd_add_member(n, &mp_glob[i]->flory_param[0][0], (MAX_CONC * MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->diffusivity_Stefan_Maxwell[0][0], (MAX_CONC * MAX_CONC),
                   MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->u_diffusivity_Stefan_Maxwell[0][0][0],
                   (MAX_CONC * MAX_CONC * MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->diffusivity_gen_fick[0][0], (MAX_CONC * MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_diffusivity_gf[0][0][0], (MAX_CONC * MAX_CONC * MAX_CONC),
                   MPI_DOUBLE);

    /*
     * Ryan's Qtensor model parameters ...
     */
    ddd_add_member(n, &mp_glob[i]->QtensorExtensionPModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->QtensorNctModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->Qtensor_Extension_P, 1, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->Qtensor_Nct, 1, MPI_DOUBLE);
    /*
     * Double derivative has lots of entries...
     */

    ddd_add_member(n, &mp_glob[i]->d_d_saturation[0][0],
                   (MAX_VARIABLE_TYPES + MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_d_cap_pres[0][0],
                   (MAX_VARIABLE_TYPES + MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC),
                   MPI_DOUBLE); /*
                                 * TFMP model parameters
                                 */
    ddd_add_member(n, &mp_glob[i]->tfmp_diff_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_diff_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->tfmp_wt_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->tfmp_wt_const, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->tfmp_density_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_density_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->tfmp_viscosity_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_viscosity_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->tfmp_rel_perm_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_rel_perm_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->tfmp_mass_lump, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->tfmp_clipping, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->tfmp_clip_strength, 1, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->tfmp_dissolution_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_dissolution_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->tfmp_drop_lattice_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_tfmp_drop_lattice_const, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->shell_tangent_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_shell_tangent_seed_vec_const, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->shell_moment_tensor_model, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->ehl_gap_model, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ehl_normal_method, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->ehl_integration_kind, 1, MPI_INT);

    /*
     * Loop over user-defined constants lists of lengths.
     *
     * The actual lists of constants will be transported via dove below.
     */

    ddd_add_member(n, &mp_glob[i]->len_u_Volume_Expansion, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_current_source, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_moment_source, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_density, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_electrical_conductivity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_permittivity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_magnetic_permeability, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_incident_wave, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_heat_capacity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_heat_source, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_mass_source, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_mesh_source, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_momentum_source, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->len_u_porosity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_porous_compressibility, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_permeability, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_rel_gas_perm, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_rel_liq_perm, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_saturation, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_cap_pres, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_tau_y, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_atexp, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_mu0, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_nexp, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_muinf, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_aexp, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_atexp, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_wlfc2, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_tau_y, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_lam, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->len_u_thixo, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->len_u_surface_tension, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_thermal_conductivity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_viscosity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_dilationalViscosity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_FlowingLiquid_viscosity, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_reaction_rate, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_solution_temperature, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_thermodynamic_potential, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_interfacial_area, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_shell_user_par, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_acoustic_impedance, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_wave_number, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_acoustic_absorption, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_refractive_index, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_light_absorption, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_extinction_index, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->thermal_conductivity_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->acoustic_impedance_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->wave_number_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->acoustic_absorption_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->refractive_index_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->light_absorption_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->extinction_index_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->viscosity_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->dilationalViscosity_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->heat_capacity_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->diffusivity_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->saturation_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->cap_pres_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LubInt_NGP, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->LubInt_PL, 1, MPI_DOUBLE);

    /*
     * Material property constants that are vectors over the concentration
     * index.
     */

    ddd_add_member(n, mp_glob[i]->diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->gam_diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->mu_diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->f_diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->g_diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->SBM_Lengths, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->NSCoeff, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->cur_diffusivity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->q_diffusivity[0][0], MAX_CONC * DIM, MPI_DOUBLE);

    ddd_add_member(n, mp_glob[i]->latent_heat_fusion, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->latent_heat_vap, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->mass_flux, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->mass_transfer_coefficient, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->reference_concn, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->species_activity, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->species_source, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->species_vol_expansion, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->vapor_pressure, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->AdvectiveScaling, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->molar_volume, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->specific_volume, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->molecular_weight, MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->charge_number, MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, mp_ptr->PhaseID, MAX_CONC, MPI_INT);
    ddd_add_member(n, &mp_ptr->Species_Var_Type, 1, MPI_INT);
    ddd_add_member(n, mp_ptr->Volumetric_Dirichlet_Cond, MAX_CONC, MPI_INT);

    ddd_add_member(n, mp_glob[i]->DiffusivityModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->LatentHeatFusionModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->LatentHeatVapModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->RefConcnModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SpecVolExpModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SpeciesSourceModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SpeciesTimeIntegration, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SpeciesOnlyDiffusion, MAX_CONC, MPI_INT);

    ddd_add_member(n, mp_glob[i]->SpeciesSecondLevelSetDiffusivity, MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, mp_glob[i]->FreeVolSolvent, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->VaporPressureModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->AdvectiveScalingModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->ExtrinsicIndependentSpeciesVar, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->MolarVolumeModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SpecificVolumeModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->MolecularWeightModel, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->ChargeNumberModel, MAX_CONC, MPI_INT);

    ddd_add_member(n, mp_glob[i]->GamDiffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->MuDiffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->GravDiffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->SBM_Type, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->NSCoeffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->FickDiffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->CurvDiffType, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->QTensorDiffType, MAX_CONC, MPI_INT);

    /*
     * Material property constants that are vectors over the porous phases
     * index.
     */

    ddd_add_member(n, mp_glob[i]->AdvectiveScaling, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->diffusivity, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->latent_heat_vap, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->latent_heat_fusion, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->vapor_pressure, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->species_vol_expansion, MAX_PMV, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->molecular_weight, MAX_PMV, MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->Porous_Mass_Lump, 1, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousDiffusivityModel, MAX_PMV, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousLatentHeatVapModel, MAX_PMV, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousLatentHeatFusionModel, MAX_PMV, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousVaporPressureModel, MAX_PMV, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousGasConstantsModel, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->PorousSinkConstantsModel, 1, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorVolExpModel, MAX_PMV, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousMolecularWeightModel, MAX_PMV, MPI_INT);

    ddd_add_member(n, mp_glob[i]->PorousShellPorosity, MAX_POR_SHELL, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->PorousShellHeight, MAX_POR_SHELL, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->PorousShellPermeability, MAX_POR_SHELL, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->PorousShellCrossPermeability, MAX_POR_SHELL, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->PorousShellRelPerm, MAX_POR_SHELL, MPI_DOUBLE);
    ddd_add_member(n, mp_glob[i]->PorousShellCapPres, MAX_POR_SHELL, MPI_DOUBLE);

    ddd_add_member(n, mp_glob[i]->PorousShellPorosityModel, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousShellHeightModel, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousShellPermeabilityModel, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousShellCrossPermeabilityModel, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousShellRelPermModel, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->PorousShellCapPresModel, MAX_POR_SHELL, MPI_INT);

    ddd_add_member(n, mp_glob[i]->por_shell_porosity_ext_field_index, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_height_ext_field_index, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_permeability_ext_field_index, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_cross_permeability_ext_field_index, MAX_POR_SHELL,
                   MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_rel_perm_ext_field_index, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_cap_pres_ext_field_index, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_cap_pres_hyst_curve_type_ext_field_index, MAX_POR_SHELL,
                   MPI_INT);
    ddd_add_member(n, mp_glob[i]->por_shell_cap_pres_hyst_num_switch_ext_field_index, MAX_POR_SHELL,
                   MPI_INT);

    /** Material property user constant lists that are vectors need
     * length variables that are vectors, too. The actual user variables
     * will be sent via the dove.
     */

    /*
     *    vectors of length to accomodate species
     */
    ddd_add_member(n, mp_glob[i]->len_u_diffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_gadiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_mdiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_cdiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_fdiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_gdiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_SBM_Lengths2, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_nscoeff, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_qdiffusivity, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_species_source, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_species_vol_expansion, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_vapor_pressure, MAX_CONC, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_reference_concn, MAX_CONC, MPI_INT);

    /*
     *    vectors of length to accomodate porous phases
     */

    ddd_add_member(n, mp_glob[i]->len_u_porous_diffusivity, MAX_PMV, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_porous_vapor_pressure, MAX_PMV, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_porous_gas_constants, 1, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->len_u_porous_sink_constants, 1, MPI_INT);

    ddd_add_member(n, mp_glob[i]->len_u_porous_vol_expansion, MAX_PMV, MPI_INT);

    ddd_add_member(n, mp_glob[i]->len_u_PorousShellPorosity, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_PorousShellHeight, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_PorousShellPermeability, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_PorousShellCrossPermeability, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_PorousShellRelPerm, MAX_POR_SHELL, MPI_INT);
    ddd_add_member(n, mp_glob[i]->len_u_PorousShellCapPres, MAX_POR_SHELL,
                   MPI_INT); /* Special material properties that are for the lubrication equation */
    ddd_add_member(n, &mp_glob[i]->len_u_heightU_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_heightL_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->heightL_function_constants_tableid, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_veloU_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_veloL_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_dcaU_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_dcaL_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_lubsource, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_lubmomsource, 1, MPI_INT);

    /* Porous shell */
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellClosedPorosity_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellClosedRadius_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellClosedHeight_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellClosedP0_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellPatm_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellPref_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellCrossKappa_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellInitPorePres_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellDiffusivity_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellRT_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_PorousShellHenry_function_constants, 1, MPI_INT);

    /* Thin film */
    ddd_add_member(n, &mp_glob[i]->len_u_FilmEvap_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_DisjPress_function_constants, 1, MPI_INT);
    ddd_add_member(n, &mp_glob[i]->len_u_DiffCoeff_function_constants, 1, MPI_INT);

    /*
     * Material properties that are fixed size arrays governed by
     * the maximum number of species and maximum number of dependent
     * variables.
     */

    ddd_add_member(n, &mp_glob[i]->d_diffusivity[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_latent_heat_fusion[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_latent_heat_vap[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_mass_flux[0][0], (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC),
                   MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_mass_transfer_coefficient[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_species_activity[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_species_vol_expansion[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_vapor_pressure[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_molecular_weight[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_charge_number[0][0],
                   (MAX_CONC) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    /*
     * Material properties that are fixed size arrays governed by
     * the maximum number of porous equations  and maximum number of dependent
     * variables.
     */

    ddd_add_member(n, &mp_glob[i]->d_PorousShellPorosity[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_PorousShellHeight[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_PorousShellPermeability[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_PorousShellCrossPermeability[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_PorousShellRelPerm[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_PorousShellCapPres[0][0],
                   (MAX_POR_SHELL) * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);

    /*
     * Material properties that are vectors with 1 dimension index.
     */

    ddd_add_member(n, &mp_glob[i]->mesh_source[0], DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->real_solid_source[0], DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->momentum_source[0], DIM, MPI_DOUBLE);

    /*
     * Material properties that are arrays with 2 dimension indeces.
     */

    ddd_add_member(n, &mp_glob[i]->perm_tensor[0][0], DIM * DIM, MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->d_perm_tensor[0][0][0],
                   DIM * DIM * (MAX_VARIABLE_TYPES + MAX_CONC), MPI_DOUBLE);
    ddd_add_member(n, &mp_glob[i]->PorousShellPermTensor[0][0][0], MAX_POR_SHELL * DIM * DIM,
                   MPI_DOUBLE);

    /*
     * Material properties that are arrays with 1 dimension and 1 vbl index.
     */

    ddd_add_member(n, &mp_glob[i]->d_mesh_source[0][0], DIM * (MAX_VARIABLE_TYPES + MAX_CONC),
                   MPI_DOUBLE);

    ddd_add_member(n, &mp_glob[i]->d_momentum_source[0][0], DIM * (MAX_VARIABLE_TYPES + MAX_CONC),
                   MPI_DOUBLE);

    /*
     * Material properties depending on variable indeces but not the
     * concentration variables !?! Are these just dinosaurs?
     */

    ddd_add_member(n, &mp_glob[i]->ReferenceModel[0], MAX_VARIABLE_TYPES, MPI_INT);

    ddd_add_member(n, &mp_glob[i]->reference[0], MAX_VARIABLE_TYPES, MPI_DOUBLE);

    /*
     * Oddball new kid on the block.
     */

    ddd_add_member(n, &mp_glob[i]->geometry_parameters[0], MAX_GEOMETRY_PARMS, MPI_DOUBLE);

    /*
     * Material properties for second level set phase
     */

    if (mp_glob[i]->mp2nd != NULL) {
      ddd_add_member(n, &mp_glob[i]->mp2nd->ViscosityModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->viscosity, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->viscositymask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->FlowingLiquidViscosityModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->FlowingLiquid_viscosity, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->FlowingLiquid_viscositymask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->DensityModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->density, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->densitymask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->HeatCapacityModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->heatcapacity, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->heatcapacitymask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->HeatSourceModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->heatsource, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->heatsourcemask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->ThermalConductivityModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->thermalconductivity, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->thermalconductivitymask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->MomentumSourceModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->momentumsource, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->momentumsourcemask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->wavenumberModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->wavenumber, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->wavenumbermask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->AcousticImpedanceModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->acousticimpedance, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->acousticimpedancemask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->AcousticAbsorptionModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->acousticabsorption, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->acousticabsorptionmask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->RefractiveIndexModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->refractiveindex, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->refractiveindexmask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->LightAbsorptionModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->lightabsorption, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->lightabsorptionmask[0], 2, MPI_INT);

      ddd_add_member(n, &mp_glob[i]->mp2nd->ExtinctionIndexModel, 1, MPI_INT);
      ddd_add_member(n, &mp_glob[i]->mp2nd->extinctionindex, 1, MPI_DOUBLE);
      ddd_add_member(n, &mp_glob[i]->mp2nd->extinctionindexmask[0], 2, MPI_INT);
      int w;
      for (w = 0; w < MAX_CONC; w++) {
        ddd_add_member(n, &mp_glob[i]->mp2nd->SpeciesSourceModel[w], 1, MPI_INT);
        ddd_add_member(n, &mp_glob[i]->mp2nd->speciessource[w], 1, MPI_DOUBLE);
        ddd_add_member(n, &mp_glob[i]->mp2nd->speciessourcemask[0][w], 1, MPI_INT);
        ddd_add_member(n, &mp_glob[i]->mp2nd->speciessourcemask[1][w], 1, MPI_INT);
        ddd_add_member(n, &mp_glob[i]->mp2nd->use_species_source_width[w], 1, MPI_INT);
        ddd_add_member(n, &mp_glob[i]->mp2nd->species_source_width[w], 1, MPI_DOUBLE);
      }
    }

    /*
     * Constitutive Relations are given for each material, too.
     */

    ddd_add_member(n, &cr_glob[i]->HeatFluxModel, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->MeshFluxModel, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->RealSolidFluxModel, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->MeshMotion, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->MassFluxModel, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->MomentumFluxModel, 1, MPI_INT);
    ddd_add_member(n, &cr_glob[i]->PorousFluxModel, 1, MPI_INT);

    /*
     * ViscoElastic models potentially apply to each material.
     */

    ddd_add_member(n, &vn_glob[i]->ConstitutiveEquation, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->ptt_type, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->wt_func, 1, MPI_DOUBLE);
    ddd_add_member(n, &vn_glob[i]->wt_funcModel, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->eps, 1, MPI_DOUBLE);
    ddd_add_member(n, &vn_glob[i]->evssModel, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->modes, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->shiftModel, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->dg_J_model, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->shockcaptureModel, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->shockcapture, 1, MPI_DOUBLE);
    /*
     * Dadblastit!!! If you add a shiny new variable argument list
     * take care to teleport it properly over to the other processors.
     *
     * Since this is now a user variable number of parameters, a length
     * variable must be declared and transported first, then, later,
     * the double vector "dg_J_model_wt" will be transported. At this
     * point, however, we have no way of knowing how much space will
     * be need for it except on processor 0 that went through mm_input.c
     */

    ddd_add_member(n, &vn_glob[i]->len_dg_J_model_wt, 1, MPI_INT);
    ddd_add_member(n, &vn_glob[i]->len_shift, 1, MPI_INT);

    for (mode = 0; mode < MAX_MODES; mode++) {
      ddd_add_member(n, &ve_glob[i][mode]->gn->ConstitutiveEquation, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->mu0, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->mu0Model, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->pos_ls_mup, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->nexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->nexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->muinf, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->muinfModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->wlfc2, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->wlfc2Model, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->lam, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->lamModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->aexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->aexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->atexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->atexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->tau_y, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->tau_yModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->fexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->fexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->epsilon, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->epsilonModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->thixo_factor, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->thixoModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->maxpack, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->maxpackModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->sus_species_no, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->gelpoint, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->gelpointModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->cureaexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->cureaexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->curebexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->curebexpModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->tgel0, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->tgel0Model, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->cure_species_no, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->gn->k1, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->k2, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->n0, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->pexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->qexp, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->diff, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->DilVisc0, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->gn->DilViscModel, 1, MPI_INT);

      ddd_add_member(n, &ve_glob[i][mode]->time_const, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->time_constModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->alpha, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->alphaModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->xi, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->xiModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->eps, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->epsModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->muJeffreys, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->muJeffreysModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->stretch_time, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->stretchModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->CCR_coefficient, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->CCR_coefficientModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->polymer_exponent, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->polymer_exponentModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->maximum_stretch_ratio, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->maximum_stretch_ratioModel, 1, MPI_INT);
      ddd_add_member(n, &ve_glob[i][mode]->extensibility, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->extensibilityModel, 1, MPI_INT);

      ddd_add_member(n, &ve_glob[i][mode]->pos_ls.time_const, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->pos_ls.alpha, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->pos_ls.xi, 1, MPI_DOUBLE);
      ddd_add_member(n, &ve_glob[i][mode]->pos_ls.eps, 1, MPI_DOUBLE);
    }

    /*
     * Generalized Newtonian structures are evidently separate entities from
     * those allocated within the context of ve_glob->
     * and might need to be transported, too.
     */

    ddd_add_member(n, &gn_glob[i]->ConstitutiveEquation, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->mu0, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->mu0Model, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->nexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->nexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->muinf, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->muinfModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->lam, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->lamModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->aexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->aexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->atexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->atexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->wlfc2, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->wlfc2Model, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->tau_y, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->tau_yModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->epsilon, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->epsilonModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->fexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->fexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->maxpack, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->maxpackModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->sus_species_no, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->gelpoint, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->gelpointModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->cureaexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->cureaexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->curebexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->curebexpModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->tgel0, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->tgel0Model, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->cure_species_no, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->k1, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->k2, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->n0, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->pexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->qexp, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->diff, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->DilVisc0, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->DilViscModel, 1, MPI_INT);
    ddd_add_member(n, &gn_glob[i]->thixo_factor, 1, MPI_DOUBLE);
    ddd_add_member(n, &gn_glob[i]->thixoModel, 1, MPI_INT);

    /*
     * Finally, the elastic constitutive models for solids and pseudosolids
     * are stored in this structure...
     */

    ddd_add_member(n, &elc_glob[i]->ConstitutiveEquation, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->lame_mu, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->lame_mu_model, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->len_u_mu, 1, MPI_INT);

    ddd_add_member(n, elc_glob[i]->d_lame_mu, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->lame_mu_tableid, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->lame_lambda, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->lame_lambda_model, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->len_u_lambda, 1, MPI_INT);

    ddd_add_member(n, elc_glob[i]->d_lame_lambda, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_glob[i]->lame_TempShift, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->lameTempShiftModel, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->len_u_lame_TempShift, 1, MPI_INT);

    ddd_add_member(n, elc_glob[i]->d_lame_TempShift, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_glob[i]->lame_TempShift_tableid, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->solid_viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->solid_viscosity_model, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->len_u_solid_viscosity, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->solid_dil_viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->solid_dil_viscosity_model, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->len_u_solid_dil_viscosity, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->bend_stiffness, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->bend_stiffness_model, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->len_u_bend_stiffness, 1, MPI_INT);

    ddd_add_member(n, elc_glob[i]->d_bend_stiffness, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_glob[i]->exten_stiffness, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->exten_stiffness_model, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->len_u_exten_stiffness, 1, MPI_INT);
    ddd_add_member(n, elc_glob[i]->d_exten_stiffness, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_glob[i]->poisson, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->poisson_model, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->len_u_poisson, 1, MPI_INT);
    ddd_add_member(n, elc_glob[i]->d_poisson, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_glob[i]->Strss_fr_sol_vol_frac, 1, MPI_DOUBLE);

    ddd_add_member(n, elc_glob[i]->v_mesh_sfs, DIM, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->v_mesh_sfs_model, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->len_u_v_mesh_sfs, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->thermal_expansion, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->thermal_expansion_model, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->len_u_thermal_expansion, 1, MPI_INT);
    ddd_add_member(n, &elc_glob[i]->thermal_expansion_tableid, 1, MPI_INT);

    ddd_add_member(n, &elc_glob[i]->solid_reference_temp, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_glob[i]->solid_reference_temp_model, 1, MPI_INT);

    /* Looks like new stuff to move to other procs... */
    ddd_add_member(n, &elc_rs_glob[i]->ConstitutiveEquation, 1, MPI_INT);
    ddd_add_member(n, &elc_rs_glob[i]->lame_mu, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->lame_mu_model, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->len_u_mu, 1, MPI_INT);

    ddd_add_member(n, elc_rs_glob[i]->d_lame_mu, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->lame_mu_tableid, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->lame_lambda, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->lame_lambda_model, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->len_u_lambda, 1, MPI_INT);

    ddd_add_member(n, elc_rs_glob[i]->d_lame_lambda, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &elc_rs_glob[i]->lame_TempShift, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->lameTempShiftModel, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->len_u_lame_TempShift, 1, MPI_INT);

    ddd_add_member(n, elc_rs_glob[i]->d_lame_TempShift, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->lame_TempShift_tableid, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->solid_viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->solid_viscosity_model, 1, MPI_INT);
    ddd_add_member(n, &elc_rs_glob[i]->len_u_solid_viscosity, 1, MPI_INT);
    ddd_add_member(n, &elc_rs_glob[i]->solid_dil_viscosity, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->solid_dil_viscosity_model, 1, MPI_INT);
    ddd_add_member(n, &elc_rs_glob[i]->len_u_solid_dil_viscosity, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->poisson, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->Strss_fr_sol_vol_frac, 1, MPI_DOUBLE);

    ddd_add_member(n, elc_rs_glob[i]->v_mesh_sfs, DIM, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->v_mesh_sfs_model, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->len_u_v_mesh_sfs, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->thermal_expansion, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->thermal_expansion_model, 1, MPI_INT);
    ddd_add_member(n, &elc_rs_glob[i]->len_u_thermal_expansion, 1, MPI_INT);

    ddd_add_member(n, &elc_rs_glob[i]->solid_reference_temp, 1, MPI_DOUBLE);
    ddd_add_member(n, &elc_rs_glob[i]->solid_reference_temp_model, 1, MPI_INT);

    /* No, really, this is it after we add the elasto-viscoplasticity struct */
    ddd_add_member(n, &evpl_glob[i]->ConstitutiveEquation, 1, MPI_INT);
    ddd_add_member(n, &evpl_glob[i]->update_flag, 1, MPI_INT);
    ddd_add_member(n, &evpl_glob[i]->plastic_mu, 1, MPI_DOUBLE);
    ddd_add_member(n, &evpl_glob[i]->plastic_mu_model, 1, MPI_INT);
    ddd_add_member(n, &evpl_glob[i]->len_u_plastic_mu, 1, MPI_INT);
    ddd_add_member(n, evpl_glob[i]->d_plastic_mu, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    ddd_add_member(n, &evpl_glob[i]->yield, 1, MPI_DOUBLE);
    ddd_add_member(n, &evpl_glob[i]->yield_model, 1, MPI_INT);
    ddd_add_member(n, &evpl_glob[i]->len_u_yield, 1, MPI_INT);
    ddd_add_member(n, evpl_glob[i]->d_yield, MAX_VARIABLE_TYPES + MAX_CONC, MPI_DOUBLE);

    /*Phil, Sam, I need help adding the global tensor beasts, since they
     *dynamically allocated */

    /*      ddd_add_member(n, Element_Blocks[i].ElemStorage[ielem].sat_curve_type, .... HOw do this
     * beast? */
  }

  /*
   * Post processing specifications, except for len_u_post_proc that was
   * transported on the first flight so that variably lengthed u_post_proc[]
   * could be moved with the main body.
   */

  ddd_add_member(n, &STREAM, 1, MPI_INT);
  ddd_add_member(n, &STREAM_NORMAL_STRESS, 1, MPI_INT);
  ddd_add_member(n, &STREAM_SHEAR_STRESS, 1, MPI_INT);
  ddd_add_member(n, &STREAM_TENSION, 1, MPI_INT);
  ddd_add_member(n, &MEAN_SHEAR, 1, MPI_INT);
  ddd_add_member(n, &PRESSURE_CONT, 1, MPI_INT);
  ddd_add_member(n, &SH_DIV_S_V_CONT, 1, MPI_INT);
  ddd_add_member(n, &SH_CURV_CONT, 1, MPI_INT);
  ddd_add_member(n, &FILL_CONT, 1, MPI_INT);
  ddd_add_member(n, &STRESS_CONT, 1, MPI_INT);
  ddd_add_member(n, &FIRST_INVAR_STRAIN, 1, MPI_INT);
  ddd_add_member(n, &SEC_INVAR_STRAIN, 1, MPI_INT);
  ddd_add_member(n, &THIRD_INVAR_STRAIN, 1, MPI_INT);

  ddd_add_member(n, &CAPILLARY_PRESSURE, 1, MPI_INT);
  ddd_add_member(n, &CONC_CONT, 1, MPI_INT);
  ddd_add_member(n, &CONDUCTION_VECTORS, 1, MPI_INT);
  ddd_add_member(n, &CURL_V, 1, MPI_INT);
  ddd_add_member(n, &DARCY_VELOCITY_GAS, 1, MPI_INT);
  ddd_add_member(n, &DARCY_VELOCITY_LIQ, 1, MPI_INT);
  ddd_add_member(n, &DIELECTROPHORETIC_FIELD, 1, MPI_INT);
  ddd_add_member(n, &DIELECTROPHORETIC_FIELD_NORM, 1, MPI_INT);
  ddd_add_member(n, &ENORMSQ_FIELD, 1, MPI_INT);
  ddd_add_member(n, &ENORMSQ_FIELD_NORM, 1, MPI_INT);
  ddd_add_member(n, &DIFFUSION_VECTORS, 1, MPI_INT);
  ddd_add_member(n, &DIFFUSION_VECTORS_POR_LIQ_GPHASE, 1, MPI_INT);
  ddd_add_member(n, &DIFFUSION_VECTORS_POR_AIR_GPHASE, 1, MPI_INT);
  ddd_add_member(n, &DIV_PVELOCITY, 1, MPI_INT);
  ddd_add_member(n, &DIV_VELOCITY, 1, MPI_INT);
  ddd_add_member(n, &DIV_TOTAL, 1, MPI_INT);
  ddd_add_member(n, &PP_Viscosity, 1, MPI_INT);
  ddd_add_member(n, &PP_FlowingLiquid_Viscosity, 1, MPI_INT);
  ddd_add_member(n, &PP_VolumeFractionGas, 1, MPI_INT);
  ddd_add_member(n, &DENSITY, 1, MPI_INT);
  ddd_add_member(n, &POLYMER_VISCOSITY, 1, MPI_INT);
  ddd_add_member(n, &POLYMER_TIME_CONST, 1, MPI_INT);
  ddd_add_member(n, &MOBILITY_PARAMETER, 1, MPI_INT);
  ddd_add_member(n, &PTT_XI, 1, MPI_INT);
  ddd_add_member(n, &PTT_EPSILON, 1, MPI_INT);
  ddd_add_member(n, &NS_RESIDUALS, 1, MPI_INT);
  ddd_add_member(n, &MM_RESIDUALS, 1, MPI_INT);
  ddd_add_member(n, &FLUXLINES, 1, MPI_INT);
  ddd_add_member(n, &ENERGY_FLUXLINES, 1, MPI_INT);
  ddd_add_member(n, &TIME_DERIVATIVES, 1, MPI_INT);
  ddd_add_member(n, &STRESS_TENSOR, 1, MPI_INT);
  ddd_add_member(n, &REAL_STRESS_TENSOR, 1, MPI_INT);
  ddd_add_member(n, &STRAIN_TENSOR, 1, MPI_INT);
  ddd_add_member(n, &EVP_DEF_GRAD_TENSOR, 1, MPI_INT);
  ddd_add_member(n, &LAGRANGE_CONVECTION, 1, MPI_INT);
  ddd_add_member(n, &POROUS_RHO_GAS_SOLVENTS, 1, MPI_INT);
  ddd_add_member(n, &POROUS_RHO_LPHASE, 1, MPI_INT);
  ddd_add_member(n, &POROUS_RHO_TOTAL_SOLVENTS, 1, MPI_INT);
  ddd_add_member(n, &POROUS_SATURATION, 1, MPI_INT);
  ddd_add_member(n, &REL_LIQ_PERM, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_VEL, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_VEL_ELSIZE, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_Q, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_Q_ELSIZE, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_P, 1, MPI_INT);
  ddd_add_member(n, &ERROR_ZZ_P_ELSIZE, 1, MPI_INT);
  ddd_add_member(n, &USER_POST, 1, MPI_INT);
  ddd_add_member(n, &ACOUSTIC_PRESSURE, 1, MPI_INT);
  ddd_add_member(n, &ACOUSTIC_PHASE_ANGLE, 1, MPI_INT);
  ddd_add_member(n, &ACOUSTIC_ENERGY_DENSITY, 1, MPI_INT);
  ddd_add_member(n, &LIGHT_INTENSITY, 1, MPI_INT);
  ddd_add_member(n, &ELECTRIC_FIELD, 1, MPI_INT);
  ddd_add_member(n, &ELECTRIC_FIELD_MAG, 1, MPI_INT);
  ddd_add_member(n, &PRINCIPAL_STRESS, 1, MPI_INT);
  ddd_add_member(n, &PRINCIPAL_REAL_STRESS, 1, MPI_INT);
  ddd_add_member(n, &LUB_HEIGHT, 1, MPI_INT);
  ddd_add_member(n, &LUB_HEIGHT_2, 1, MPI_INT);
  ddd_add_member(n, &LUB_VELO_UPPER, 1, MPI_INT);
  ddd_add_member(n, &LUB_VELO_LOWER, 1, MPI_INT);
  ddd_add_member(n, &LUB_VELO_FIELD, 1, MPI_INT);
  ddd_add_member(n, &LUB_VELO_FIELD_2, 1, MPI_INT);
  ddd_add_member(n, &DISJ_PRESS, 1, MPI_INT);
  ddd_add_member(n, &SH_SAT_OPEN, 1, MPI_INT);
  ddd_add_member(n, &SH_SAT_OPEN_2, 1, MPI_INT);
  ddd_add_member(n, &SH_CAP_PRES, 1, MPI_INT);
  ddd_add_member(n, &SH_PORE_FLUX, 1, MPI_INT);
  ddd_add_member(n, &SH_STRESS_TENSOR, 1, MPI_INT);
  ddd_add_member(n, &SH_TANG, 1, MPI_INT);
  ddd_add_member(n, &PP_LAME_MU, 1, MPI_INT);
  ddd_add_member(n, &PP_LAME_LAMBDA, 1, MPI_INT);
  ddd_add_member(n, &VON_MISES_STRAIN, 1, MPI_INT);
  ddd_add_member(n, &VON_MISES_STRESS, 1, MPI_INT);
  ddd_add_member(n, &UNTRACKED_SPEC, 1, MPI_INT);
  ddd_add_member(n, &CONF_MAP, 1, MPI_INT);
  ddd_add_member(n, &VELO_SPEED, 1, MPI_INT);
  ddd_add_member(n, &GIES_CRIT, 1, MPI_INT);
  ddd_add_member(n, &J_FLUX, 1, MPI_INT);
  ddd_add_member(n, &EIG, 1, MPI_INT);
  ddd_add_member(n, &HEAVISIDE, 1, MPI_INT);
  ddd_add_member(n, &RHO_DOT, 1, MPI_INT);
  ddd_add_member(n, &MOMENT_SOURCES, 1, MPI_INT);
  ddd_add_member(n, &YZBETA, 1, MPI_INT);
  ddd_add_member(n, &EIG1, 1, MPI_INT);
  ddd_add_member(n, &EIG2, 1, MPI_INT);
  ddd_add_member(n, &EIG3, 1, MPI_INT);
  ddd_add_member(n, &GRAD_SH, 1, MPI_INT);
  ddd_add_member(n, &GRAD_Y, 1, MPI_INT);
  ddd_add_member(n, &HELICITY, 1, MPI_INT);
  ddd_add_member(n, &EXTERNAL_POST, 1, MPI_INT);
  ddd_add_member(n, &SURFACE_VECTORS, 1, MPI_INT);
  ddd_add_member(n, &SHELL_NORMALS, 1, MPI_INT);
  ddd_add_member(n, &len_u_post_proc, 1, MPI_INT);
  ddd_add_member(n, &LAMB_VECTOR, 1, MPI_INT);
  ddd_add_member(n, &Q_FCN, 1, MPI_INT);
  ddd_add_member(n, &POYNTING_VECTORS, 1, MPI_INT);
  ddd_add_member(n, &SARAMITO_YIELD, 1, MPI_INT);
  ddd_add_member(n, &STRESS_NORM, 1, MPI_INT);
  ddd_add_member(n, &SPECIES_SOURCES, 1, MPI_INT);
  ddd_add_member(n, &VISCOUS_STRESS, 1, MPI_INT);
  ddd_add_member(n, &VISCOUS_STRESS_NORM, 1, MPI_INT);
  ddd_add_member(n, &VISCOUS_VON_MISES_STRESS, 1, MPI_INT);
  ddd_add_member(n, &EM_CONTOURS, 1, MPI_INT);
  ddd_add_member(n, &TOTAL_EM_CONTOURS, 1, MPI_INT);
  ddd_add_member(n, &SCATTERED_EM_CONTOURS, 1, MPI_INT);
  ddd_add_member(n, &len_u_post_proc, 1, MPI_INT);
  ddd_add_member(n, &PSPG_PP, 1, MPI_INT);
  ddd_add_member(n, &ORIENTATION_VECTORS, 1, MPI_INT);
  ddd_add_member(n, &FIRST_STRAINRATE_INVAR, 1, MPI_INT);
  ddd_add_member(n, &SEC_STRAINRATE_INVAR, 1, MPI_INT);
  ddd_add_member(n, &THIRD_STRAINRATE_INVAR, 1, MPI_INT);
  ddd_add_member(n, &WALL_DISTANCE, 1, MPI_INT);
  if (len_u_post_proc > 0) {
    ddd_add_member(n, u_post_proc, len_u_post_proc, MPI_DOUBLE);
  }

  if (nn_post_fluxes > 0) {
    for (i = 0; i < nn_post_fluxes; i++) {
      ddd_add_member(n, &(pp_fluxes[i]->flux_type), 1, MPI_INT);
      ddd_add_member(n, pp_fluxes[i]->flux_type_name, MAX_DOFNAME, MPI_CHAR);
      ddd_add_member(n, &(pp_fluxes[i]->species_number), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes[i]->ss_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes[i]->blk_id), 1, MPI_INT);
      ddd_add_member(n, pp_fluxes[i]->flux_filenm, MAX_FNL, MPI_CHAR);
      ddd_add_member(n, &(pp_fluxes[i]->profile_flag), 1, MPI_INT);
    }
  }

  if (nn_post_fluxes_sens > 0) {
    for (i = 0; i < nn_post_fluxes_sens; i++) {
      ddd_add_member(n, &(pp_fluxes_sens[i]->flux_type), 1, MPI_INT);
      ddd_add_member(n, pp_fluxes_sens[i]->flux_type_name, MAX_DOFNAME, MPI_CHAR);
      ddd_add_member(n, &(pp_fluxes_sens[i]->species_number), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->ss_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->blk_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->sens_type), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->sens_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->sens_flt), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->sens_flt2), 1, MPI_INT);
      ddd_add_member(n, &(pp_fluxes_sens[i]->vector_id), 1, MPI_INT);
      ddd_add_member(n, pp_fluxes_sens[i]->flux_filenm, MAX_FNL, MPI_CHAR);
      ddd_add_member(n, &(pp_fluxes_sens[i]->profile_flag), 1, MPI_INT);
    }
  }

  if (nn_post_data > 0) {
    for (i = 0; i < nn_post_data; i++) {
      ddd_add_member(n, &(pp_data[i]->data_type), 1, MPI_INT);
      ddd_add_member(n, pp_data[i]->data_type_name, MAX_DOFNAME, MPI_CHAR);
      ddd_add_member(n, &(pp_data[i]->species_number), 1, MPI_INT);
      ddd_add_member(n, &(pp_data[i]->ns_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_data[i]->mat_num), 1, MPI_INT);
      ddd_add_member(n, &(pp_data[i]->elem_blk_id), 1, MPI_INT);
      ddd_add_member(n, pp_data[i]->data_filenm, MAX_FNL, MPI_CHAR);
      ddd_add_member(n, pp_data[i]->format_flag, 8, MPI_CHAR);
      ddd_add_member(n, &(pp_data[i]->first_time), 1, MPI_INT);
    }
  }

  if (nn_post_data_sens > 0) {
    for (i = 0; i < nn_post_data_sens; i++) {
      ddd_add_member(n, &(pp_data_sens[i]->data_type), 1, MPI_INT);
      ddd_add_member(n, pp_data_sens[i]->data_type_name, MAX_DOFNAME, MPI_CHAR);
      ddd_add_member(n, &(pp_data_sens[i]->species_number), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->ns_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->mat_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->sens_type), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->sens_id), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->sens_flt), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->sens_flt2), 1, MPI_INT);
      ddd_add_member(n, &(pp_data_sens[i]->vector_id), 1, MPI_INT);
      ddd_add_member(n, pp_data_sens[i]->data_filenm, MAX_FNL, MPI_CHAR);
    }
  }

  if (nn_volume > 0) {
    for (i = 0; i < nn_volume; i++) {
      ddd_add_member(n, &(pp_volume[i]->volume_type), 1, MPI_INT);
      ddd_add_member(n, pp_volume[i]->volume_name, MAX_DOFNAME, MPI_CHAR);
      ddd_add_member(n, &(pp_volume[i]->species_no), 1, MPI_INT);
      ddd_add_member(n, &(pp_volume[i]->blk_id), 1, MPI_INT);
      ddd_add_member(n, pp_volume[i]->volume_fname, MAX_FNL, MPI_CHAR);
      ddd_add_member(n, &(pp_volume[i]->num_params), 1, MPI_INT);
    }
  }

  if (nn_global > 0) {
    for (i = 0; i < nn_global; i++) {
      ddd_add_member(n, &(pp_global[i]->type), 1, MPI_INT);
      ddd_add_member(n, pp_global[i]->filenm, MAX_FNL, MPI_CHAR);
    }
  }

  if (nn_average > 0) {
    for (i = 0; i < nn_average; i++) {
      ddd_add_member(n, &(pp_average[i]->type), 1, MPI_INT);
      ddd_add_member(n, &(pp_average[i]->species_index), 1, MPI_INT);
      ddd_add_member(n, &(pp_average[i]->index), 1, MPI_INT);
      ddd_add_member(n, &(pp_average[i]->index_post), 1, MPI_INT);
      ddd_add_member(n, &(pp_average[i]->non_variable_type), 1, MPI_INT);
      ddd_add_member(n, pp_average[i]->type_name, MAX_VAR_NAME_LNGTH, MPI_CHAR);
    }
  }

/*
  if ( nn_error_metrics > 0 )
    {
      for ( i=0; i < nn_error_metrics; i++ )
        {
          ddd_add_member(n, &pp_error_data[i], 6, MPI_DOUBLE);
        }
    }
*/
#endif

  ddd_set_commit(n);

  return;
}

/*
 * ark_landing() -- polish off some tasks on processors other
 *                  than processor zero, after the ark.
 */

void ark_landing(void) {
  int i;
  int index;
  int j;
  struct Elastic_Constitutive *e;
  struct Material_Properties *m;
  struct Viscoelastic_Nonmodal *v;

  if (ProcID == 0)
    return;

  // turbulence information
  if (upd->turbulent_info->num_node_sets > 0) {
    upd->turbulent_info->node_set_ids =
        (int *)malloc(upd->turbulent_info->num_node_sets * sizeof(int));
  }
  if (upd->turbulent_info->num_side_sets > 0) {
    upd->turbulent_info->side_set_ids =
        (int *)malloc(upd->turbulent_info->num_side_sets * sizeof(int));
  }

  /* Particle variable-lengthed output variable lists. */
  for (i = 0; i < Particle_Number_Sample_Types; i++) {
    if (Particle_Number_Output_Variables[i])
      Particle_Output_Variables[i] = (particle_variable_s *)array_alloc(
          1, Particle_Number_Output_Variables[i], sizeof(particle_variable_s));
    else
      Particle_Output_Variables[i] = NULL;
  }

  /*
   * Set up some pointers based on the indeces received from Processor 0
   * so they make sense here...
   *
   * Example: the pointers to BC Descriptions, static. Use the index.
   */

  for (i = 0; i < Num_BC; i++) {

    /*
     * Setup pointers into the static BC Descriptions for most BCs.
     */

    if (BC_Types[i].index_dad == -1) {
      BC_Types[i].desc = &BC_Desc[BC_Types[i].BC_Desc_index];
    } else {
      BC_Types[i].desc = NULL;
    }

    /*
     * BCs requiring a new dynamic BC Description obtain pointers
     * to old names (just like what happens in alloc_BC_Desc().)
     */

    if (BC_Types[i].index_dad != -1) {
      BC_Types[i].desc = new_BC_Desc[BC_Types[i].index_dad];
      index = BC_Types[i].BC_Desc_index;
      BC_Types[i].desc->name1 = BC_Desc[index].name1;
      BC_Types[i].desc->name2 = BC_Desc[index].name2;
    }

    /*
     * Finally, any BC's that have nontrivial numbers of user-defined
     * double constants must have space allocated on processors other
     * than ProcID 0.
     */

    dalloc(BC_Types[i].len_u_BC, BC_Types[i].u_BC);

    /*
     * Handle BC Tables
     */

    if (BC_Types[i].table_index != -1) {
      BC_Types[i].table = BC_Tables[BC_Types[i].table_index];
    } else {
      BC_Types[i].table = NULL;
    }
  }

  /*
   * Allocate arrays for the BC_Table data,
   *  note that there is no need to assign
   *  the axis (t) or (f) name on all procs
   */

  for (i = 0; i < num_BC_Tables; i++) {
    dalloc(BC_Tables[i]->tablelength, BC_Tables[i]->t);
    dalloc(BC_Tables[i]->tablelength, BC_Tables[i]->t2);
    dalloc(BC_Tables[i]->tablelength, BC_Tables[i]->f);
  }

  /*
   * Allocate space for any defined ext_Tables
   */

  for (i = 0; i < num_ext_Tables; i++) {
    dalloc(ext_Tables[i]->tablelength, ext_Tables[i]->t);
    dalloc(ext_Tables[i]->tablelength, ext_Tables[i]->t2);
    dalloc(ext_Tables[i]->tablelength, ext_Tables[i]->t3);
    dalloc(ext_Tables[i]->tablelength, ext_Tables[i]->f);
  }

  /*
   * Reconstruct pointers to BC_descriptions embedded in the freshly arrived
   * Rotation structures.
   *
   * Unlike the BC_Types, the ROT_Types refer solely to existing BC_Desc's
   * and, not yet, any dynamically allocated BC_Descs. If they ever do, then
   * a similar procedure for transporting references to dynamically allocated
   * memory across processors would need to be carried out (i.e., add a
   * "int index_dad[CDIM]" to the Rotation_Specs structure...)
   */

  for (i = 0; i < Num_ROT; i++) {
    for (j = 0; j < CDIM; j++) {
      index = ROT_Types[i].BC_desc_index[j];
      if (index > -1) {
        ROT_Types[i].BC_desc[j] = &BC_Desc[index];
      } else {
        ROT_Types[i].BC_desc[j] = NULL;
      }
    }
  }

  /*
   * Allocate space for any defined MP_Tables
   */

  for (i = 0; i < num_MP_Tables; i++) {
    dalloc(MP_Tables[i]->tablelength, MP_Tables[i]->t);
    dalloc(MP_Tables[i]->tablelength, MP_Tables[i]->t2);
    dalloc(MP_Tables[i]->tablelength, MP_Tables[i]->f);
  }

  /*
   * Conditionally allocate space for any user-defined lists of double
   * constants, as long as the length is nontrivial. Likewise, the dove
   * will register the types of each nontrivial list with MPI for transport.
   */

  for (i = 0; i < upd->Num_Mat; i++) {

    /*
     * Material properties...
     */

    m = mp_glob[i];

    /*
     *  Allocate storage in the material structure based upon
     *  the number of element blocks in the material.
     */
    m->Matrl_Elem_Blk_Ids = alloc_int_1(m->Num_Matrl_Elem_Blk, -1);

    /*
     *  Allocate space based on the number of species in the
     *  material
     */
    m->Species_Names = alloc_VecFixedStrings(m->Num_Species, sizeof(CK_NAME_STR));

    /*
     *  Allocate space based on the number of porous media phases
     *  in the material
     */
    m->Porous_Names = alloc_VecFixedStrings(m->Num_Porous_Eqn, sizeof(CK_NAME_STR));

    dalloc(m->len_u_Volume_Expansion, m->u_Volume_Expansion);

    dalloc(m->len_u_current_source, m->u_current_source);

    dalloc(m->len_u_moment_source, m->u_moment_source);

    dalloc(m->len_u_density, m->u_density);

    dalloc(m->len_u_electrical_conductivity, m->u_electrical_conductivity);

    dalloc(m->len_u_permittivity, m->u_permittivity);

    dalloc(m->len_u_magnetic_permeability, m->u_magnetic_permeability);

    dalloc(m->len_u_incident_wave, m->u_incident_wave);

    dalloc(m->len_u_elect_surf_diffusivity, m->u_elect_surf_diffusivity);

    dalloc(m->len_u_heat_capacity, m->u_heat_capacity);

    dalloc(m->len_u_heat_source, m->u_heat_source);

    dalloc(m->len_u_mass_source, m->u_mass_source);

    dalloc(m->len_u_mesh_source, m->u_mesh_source);

    dalloc(m->len_u_momentum_source, m->u_momentum_source);

    dalloc(m->len_u_porosity, m->u_porosity);

    dalloc(m->len_u_porous_compressibility, m->u_porous_compressibility);

    dalloc(m->len_u_permeability, m->u_permeability);

    dalloc(m->len_u_rel_gas_perm, m->u_rel_gas_perm);

    dalloc(m->len_u_rel_liq_perm, m->u_rel_liq_perm);

    dalloc(m->len_u_saturation, m->u_saturation);

    dalloc(m->len_u_cap_pres, m->u_cap_pres);
    dalloc(m->len_u_surface_tension, m->u_surface_tension);

    dalloc(m->len_u_thermal_conductivity, m->u_thermal_conductivity);

    dalloc(m->len_u_viscosity, m->u_viscosity);

    dalloc(m->len_u_dilationalViscosity, m->u_dilationalViscosity);

    dalloc(gn_glob[i]->len_u_mu0, gn_glob[i]->u_mu0);

    dalloc(gn_glob[i]->len_u_nexp, gn_glob[i]->u_nexp);

    dalloc(gn_glob[i]->len_u_muinf, gn_glob[i]->u_muinf);

    dalloc(gn_glob[i]->len_u_aexp, gn_glob[i]->u_aexp);

    dalloc(gn_glob[i]->len_u_atexp, gn_glob[i]->u_atexp);

    dalloc(gn_glob[i]->len_u_wlfc2, gn_glob[i]->u_wlfc2);

    dalloc(gn_glob[i]->len_u_lam, gn_glob[i]->u_lam);

    dalloc(gn_glob[i]->len_u_tau_y, gn_glob[i]->u_tau_y);

    dalloc(gn_glob[i]->len_u_thixo, gn_glob[i]->u_thixo_factor);

    dalloc(m->len_u_FlowingLiquid_viscosity, m->u_FlowingLiquid_viscosity);

    dalloc(m->len_u_reaction_rate, m->u_reaction_rate);

    dalloc(m->len_u_solution_temperature, m->u_solution_temperature);

    dalloc(m->len_u_thermodynamic_potential, m->u_thermodynamic_potential);

    dalloc(m->len_u_interfacial_area, m->u_interfacial_area);

    dalloc(m->len_u_shell_user_par, m->u_shell_user_par);

    dalloc(m->len_u_acoustic_impedance, m->u_acoustic_impedance);

    dalloc(m->len_u_wave_number, m->u_wave_number);

    dalloc(m->len_u_acoustic_absorption, m->u_acoustic_absorption);

    dalloc(m->len_u_refractive_index, m->u_refractive_index);

    dalloc(m->len_u_light_absorption, m->u_light_absorption);

    dalloc(m->len_u_extinction_index, m->u_extinction_index);

    dalloc(m->len_tfmp_density_const, m->tfmp_density_const);
    dalloc(m->len_tfmp_viscosity_const, m->tfmp_viscosity_const);
    dalloc(m->len_tfmp_diff_const, m->tfmp_diff_const);

    dalloc(m->len_tfmp_rel_perm_const, m->tfmp_rel_perm_const);
    dalloc(m->len_tfmp_dissolution_const, m->tfmp_dissolution_const);
    dalloc(m->len_tfmp_drop_lattice_const, m->tfmp_drop_lattice_const);

    dalloc(m->len_shell_tangent_seed_vec_const, m->shell_tangent_seed_vec_const);

    /*
     * User defined material property lists for each species...
     *     HKM -> Changed this to number of species, not
     *            number of species equations.
     */

    for (j = 0; j < pd_glob[i]->Num_Species; j++) {

      dalloc(m->len_u_diffusivity[j], m->u_diffusivity[j]);

      dalloc(m->len_u_gadiffusivity[j], m->u_gadiffusivity[j]);

      dalloc(m->len_u_cdiffusivity[j], m->u_cdiffusivity[j]);

      dalloc(m->len_u_mdiffusivity[j], m->u_mdiffusivity[j]);

      dalloc(m->len_u_fdiffusivity[j], m->u_fdiffusivity[j]);

      dalloc(m->len_u_gdiffusivity[j], m->u_gdiffusivity[j]);

      dalloc(m->len_SBM_Lengths2[j], m->SBM_Lengths2[j]);

      dalloc(m->len_u_nscoeff[j], m->u_nscoeff[j]);

      dalloc(m->len_u_qdiffusivity[j], m->u_qdiffusivity[j]);

      dalloc(m->len_u_species_source[j], m->u_species_source[j]);

      dalloc(m->len_u_species_vol_expansion[j], m->u_species_vol_expansion[j]);

      dalloc(m->len_u_vapor_pressure[j], m->u_vapor_pressure[j]);

      dalloc(m->len_u_reference_concn[j], m->u_reference_concn[j]);
    }

    /*
     * User defined material property lists for each porous phase...
     */

    for (j = 0; j < pd_glob[i]->Num_Porous_Eqn; j++) {
      dalloc(m->len_u_porous_diffusivity[j], m->u_porous_diffusivity[j]);
      dalloc(m->len_u_porous_vol_expansion[j], m->u_porous_vol_expansion[j]);

      dalloc(m->len_u_porous_vapor_pressure[j], m->u_porous_vapor_pressure[j]);
    }
    for (j = 0; j < pd_glob[i]->Num_Porous_Shell_Eqn; j++) {
      dalloc(m->len_u_PorousShellPorosity[j], m->u_PorousShellPorosity[j]);
      dalloc(m->len_u_PorousShellHeight[j], m->u_PorousShellHeight[j]);
      dalloc(m->len_u_PorousShellPermeability[j], m->u_PorousShellPermeability[j]);
      dalloc(m->len_u_PorousShellCrossPermeability[j], m->u_PorousShellCrossPermeability[j]);
      dalloc(m->len_u_PorousShellRelPerm[j], m->u_PorousShellRelPerm[j]);
      dalloc(m->len_u_PorousShellCapPres[j], m->u_PorousShellCapPres[j]);
    }
    dalloc(m->len_u_porous_gas_constants, m->u_porous_gas_constants);

    dalloc(m->len_u_porous_sink_constants, m->u_porous_sink_constants);

    /*
     * special usage models for lubrication
     */
    dalloc(m->len_u_heightU_function_constants, m->u_heightU_function_constants);
    dalloc(m->len_u_heightL_function_constants, m->u_heightL_function_constants);
    dalloc(m->len_u_veloU_function_constants, m->u_veloU_function_constants);
    dalloc(m->len_u_veloL_function_constants, m->u_veloL_function_constants);
    dalloc(m->len_u_dcaU_function_constants, m->u_dcaU_function_constants);
    dalloc(m->len_u_dcaL_function_constants, m->u_dcaL_function_constants);
    dalloc(m->len_lubsource, m->u_lubsource_function_constants);
    /*
     * special usage models for porous shell
     */
    dalloc(m->len_u_PorousShellClosedPorosity_function_constants,
           m->u_PorousShellClosedPorosity_function_constants);
    dalloc(m->len_u_PorousShellClosedRadius_function_constants,
           m->u_PorousShellClosedRadius_function_constants);
    dalloc(m->len_u_PorousShellClosedHeight_function_constants,
           m->u_PorousShellClosedHeight_function_constants);
    dalloc(m->len_u_PorousShellClosedP0_function_constants,
           m->u_PorousShellClosedP0_function_constants);
    dalloc(m->len_u_PorousShellPatm_function_constants, m->u_PorousShellPatm_function_constants);
    dalloc(m->len_u_PorousShellPref_function_constants, m->u_PorousShellPref_function_constants);
    dalloc(m->len_u_PorousShellCrossKappa_function_constants,
           m->u_PorousShellCrossKappa_function_constants);
    dalloc(m->len_u_PorousShellInitPorePres_function_constants,
           m->u_PorousShellInitPorePres_function_constants);
    dalloc(m->len_u_PorousShellDiffusivity_function_constants,
           m->u_PorousShellDiffusivity_function_constants);
    dalloc(m->len_u_PorousShellRT_function_constants, m->u_PorousShellRT_function_constants);
    dalloc(m->len_u_PorousShellHenry_function_constants, m->u_PorousShellHenry_function_constants);
    /*
     * special usage models for thin film
     */
    dalloc(m->len_u_FilmEvap_function_constants, m->u_FilmEvap_function_constants);
    dalloc(m->len_u_DisjPress_function_constants, m->u_DisjPress_function_constants);
    dalloc(m->len_u_DiffCoeff_function_constants, m->u_DiffCoeff_function_constants);
    /*
     * Elastic Constitutive properties...
     */
    e = elc_glob[i];

    dalloc(e->len_u_mu, e->u_mu);

    dalloc(e->len_u_lambda, e->u_lambda);

    dalloc(e->len_u_lame_TempShift, e->u_lame_TempShift);

    dalloc(e->len_u_v_mesh_sfs, e->u_v_mesh_sfs);

    dalloc(e->len_u_thermal_expansion, e->u_thermal_expansion);

    dalloc(e->len_u_solid_viscosity, e->u_solid_viscosity);

    dalloc(e->len_u_solid_dil_viscosity, e->u_solid_dil_viscosity);

    e = elc_rs_glob[i];

    dalloc(e->len_u_mu, e->u_mu);

    dalloc(e->len_u_lambda, e->u_lambda);

    dalloc(e->len_u_lame_TempShift, e->u_lame_TempShift);

    dalloc(e->len_u_v_mesh_sfs, e->u_v_mesh_sfs);

    dalloc(e->len_u_thermal_expansion, e->u_thermal_expansion);

    dalloc(e->len_u_solid_viscosity, e->u_solid_viscosity);

    dalloc(e->len_u_solid_dil_viscosity, e->u_solid_dil_viscosity);

    dalloc(evpl_glob[i]->len_u_plastic_mu, evpl_glob[i]->u_plastic_mu);

    dalloc(evpl_glob[i]->len_u_yield, evpl_glob[i]->u_yield);

    /*
     * Viscoelastic nonmodal properties of variable length...
     */

    v = vn_glob[i];

    dalloc(v->len_dg_J_model_wt, v->dg_J_model_wt);
    dalloc(v->len_shift, v->shift);

    for (j = 0; j < Num_Var_Init_Mat[i]; j++) {
      dalloc(Var_init_mat[i][j].len_u_pars, Var_init_mat[i][j].u_pars);
    }
  }

  dalloc(len_u_post_proc, u_post_proc);

  /*   Augmenting condition data tables    */

  for (i = 0; i < nAC; i++) {
    dalloc(augc[i].len_AC, augc[i].DataFlt);
    if (augc[i].Aprepro_lib_string_len > 0) {
      augc[i].Aprepro_lib_string = calloc(augc[i].Aprepro_lib_string_len + 1, sizeof(char));
    }
  }

  for (i = 0; i < nn_volume; i++) {
    dalloc(pp_volume[i]->num_params, pp_volume[i]->params);
  }

  return;
}

/*
 * Dove is the 3rd & last delivery, after the 1st raven and the 2nd ark.
 *
 * Dove carries things of lengths that could not be known until now.
 *
 * Examples: variable numbers of user parameters for boundary conditions
 * and thermophysical properties are registered here for transport.
 */

void noahs_dove(void) {
#ifdef PARALLEL
  int i, j, k;
  DDD n;
  struct Elastic_Constitutive *e;
  struct Material_Properties *m;
  struct Viscoelastic_Nonmodal *v;

  Noahs_Dove = ddd_alloc();
  n = Noahs_Dove;

  for (i = 0; i < upd->Num_Mat; i++) {

    /*
     * Material properties...
     */
    m = mp_glob[i];

    if (m->Num_Matrl_Elem_Blk > 0) {
      if (m->Matrl_Elem_Blk_Ids == NULL) {
        printf("P_%d: ERROR: NULL address but nonzero length, %s line %d!\n", ProcID, __FILE__,
               __LINE__);
        fflush(stdout);
        GOMA_EH(GOMA_ERROR, "noahs_dove fatal error");
      }
      ddd_add_member(n, m->Matrl_Elem_Blk_Ids, m->Num_Matrl_Elem_Blk, MPI_INT);
    }

    crdv(m->len_u_Volume_Expansion, m->u_Volume_Expansion);

    crdv(m->len_u_current_source, m->u_current_source);

    crdv(m->len_u_moment_source, m->u_moment_source);

    crdv(m->len_u_density, m->u_density);

    crdv(m->len_u_electrical_conductivity, m->u_electrical_conductivity);

    crdv(m->len_u_permittivity, m->u_permittivity);

    crdv(m->len_u_magnetic_permeability, m->u_magnetic_permeability);

    crdv(m->len_u_incident_wave, m->u_incident_wave);

    crdv(m->len_u_elect_surf_diffusivity, m->u_elect_surf_diffusivity);

    crdv(m->len_u_heat_capacity, m->u_heat_capacity);

    crdv(m->len_u_heat_source, m->u_heat_source);

    crdv(m->len_u_mass_source, m->u_mass_source);

    crdv(m->len_u_mesh_source, m->u_mesh_source);

    crdv(m->len_u_momentum_source, m->u_momentum_source);

    crdv(m->len_u_porosity, m->u_porosity);

    crdv(m->len_u_permeability, m->u_permeability);

    crdv(m->len_u_rel_gas_perm, m->u_rel_gas_perm);

    crdv(m->len_u_rel_liq_perm, m->u_rel_liq_perm);

    crdv(m->len_u_cap_pres, m->u_cap_pres);

    crdv(m->len_u_saturation, m->u_saturation);

    crdv(m->len_u_porous_sink_constants, m->u_porous_sink_constants);

    crdv(m->len_u_surface_tension, m->u_surface_tension);

    crdv(m->len_u_thermal_conductivity, m->u_thermal_conductivity);

    crdv(m->len_u_viscosity, m->u_viscosity);

    crdv(m->len_u_dilationalViscosity, m->u_dilationalViscosity);

    crdv(gn_glob[i]->len_u_mu0, gn_glob[i]->u_mu0);

    crdv(gn_glob[i]->len_u_nexp, gn_glob[i]->u_nexp);

    crdv(gn_glob[i]->len_u_muinf, gn_glob[i]->u_muinf);

    crdv(gn_glob[i]->len_u_aexp, gn_glob[i]->u_aexp);

    crdv(gn_glob[i]->len_u_atexp, gn_glob[i]->u_atexp);

    crdv(gn_glob[i]->len_u_wlfc2, gn_glob[i]->u_wlfc2);

    crdv(gn_glob[i]->len_u_lam, gn_glob[i]->u_lam);

    crdv(gn_glob[i]->len_u_tau_y, gn_glob[i]->u_tau_y);

    crdv(gn_glob[i]->len_u_thixo, gn_glob[i]->u_thixo_factor);

    crdv(m->len_u_FlowingLiquid_viscosity, m->u_FlowingLiquid_viscosity);

    crdv(m->len_u_reaction_rate, m->u_reaction_rate);

    crdv(m->len_u_solution_temperature, m->u_solution_temperature);

    crdv(m->len_u_thermodynamic_potential, m->u_thermodynamic_potential);

    crdv(m->len_u_interfacial_area, m->u_interfacial_area);

    crdv(m->len_u_shell_user_par, m->u_shell_user_par);

    crdv(m->len_u_acoustic_impedance, m->u_acoustic_impedance);

    crdv(m->len_u_wave_number, m->u_wave_number);

    crdv(m->len_u_acoustic_absorption, m->u_acoustic_absorption);

    crdv(m->len_u_refractive_index, m->u_refractive_index);

    crdv(m->len_u_light_absorption, m->u_light_absorption);

    crdv(m->len_u_extinction_index, m->u_extinction_index);

    crdv(m->len_tfmp_density_const, m->tfmp_density_const);
    crdv(m->len_tfmp_viscosity_const, m->tfmp_viscosity_const);
    crdv(m->len_tfmp_diff_const, m->tfmp_diff_const);

    crdv(m->len_tfmp_rel_perm_const, m->tfmp_rel_perm_const);

    crdv(m->len_tfmp_drop_lattice_const, m->tfmp_drop_lattice_const);

    crdv(m->len_shell_tangent_seed_vec_const, m->shell_tangent_seed_vec_const);

    crdv(m->LubInt_NGP, m->Lub_gpts);
    crdv(m->LubInt_NGP, m->Lub_wts);

    /*
     *  Add species names
     */
    for (k = 0; k < m->Num_Species; k++) {
      ddd_add_member(n, m->Species_Names[k], sizeof(CK_NAME_STR), MPI_CHAR);
    }

    /*
     *  Add porous phase names
     */
    for (k = 0; k < m->Num_Porous_Eqn; k++) {
      ddd_add_member(n, m->Porous_Names[k], sizeof(CK_NAME_STR), MPI_CHAR);
    }

    /* Particle variable-lengthed output variable lists. */
    for (j = 0; j < Particle_Number_Sample_Types; j++)
      for (k = 0; k < Particle_Number_Output_Variables[j]; k++)
        ddd_add_member(n, Particle_Output_Variables[j][k], MAX_PARTICLE_OUTPUT_VARIABLE_LENGTH,
                       MPI_CHAR);

    /*
     * Species constants...
     */

    for (j = 0; j < pd_glob[i]->Num_Species; j++) {
      crdv(m->len_u_diffusivity[j], m->u_diffusivity[j]);
      crdv(m->len_u_gadiffusivity[j], m->u_gadiffusivity[j]);
      crdv(m->len_u_mdiffusivity[j], m->u_mdiffusivity[j]);
      crdv(m->len_u_cdiffusivity[j], m->u_cdiffusivity[j]);
      crdv(m->len_u_fdiffusivity[j], m->u_fdiffusivity[j]);
      crdv(m->len_u_gdiffusivity[j], m->u_gdiffusivity[j]);
      crdv(m->len_SBM_Lengths2[j], m->SBM_Lengths2[j]);
      crdv(m->len_u_nscoeff[j], m->u_nscoeff[j]);
      crdv(m->len_u_qdiffusivity[j], m->u_qdiffusivity[j]);
      crdv(m->len_u_species_source[j], m->u_species_source[j]);
      crdv(m->len_u_species_vol_expansion[j], m->u_species_vol_expansion[j]);
      crdv(m->len_u_vapor_pressure[j], m->u_vapor_pressure[j]);
      crdv(m->len_u_reference_concn[j], m->u_reference_concn[j]);
    }

    /*
     * Porous media user constants...
     */

    for (j = 0; j < pd_glob[i]->Num_Porous_Eqn; j++) {
      crdv(m->len_u_porous_diffusivity[j], m->u_porous_diffusivity[j]);
      crdv(m->len_u_porous_vol_expansion[j], m->u_porous_vol_expansion[j]);
      crdv(m->len_u_porous_vapor_pressure[j], m->u_porous_vapor_pressure[j]);
    }
    crdv(m->len_u_porous_gas_constants, m->u_porous_gas_constants);

    e = elc_glob[i];

    for (j = 0; j < pd_glob[i]->Num_Porous_Shell_Eqn; j++) {
      crdv(m->len_u_PorousShellPorosity[j], m->u_PorousShellPorosity[j]);
      crdv(m->len_u_PorousShellHeight[j], m->u_PorousShellHeight[j]);
      crdv(m->len_u_PorousShellPermeability[j], m->u_PorousShellPermeability[j]);
      crdv(m->len_u_PorousShellCrossPermeability[j], m->u_PorousShellCrossPermeability[j]);
      crdv(m->len_u_PorousShellRelPerm[j], m->u_PorousShellRelPerm[j]);
      crdv(m->len_u_PorousShellCapPres[j], m->u_PorousShellCapPres[j]);
    } /*
       * some more specialized constants (Lubrication)
       */
    crdv(m->len_u_heightU_function_constants, m->u_heightU_function_constants);
    crdv(m->len_u_heightL_function_constants, m->u_heightL_function_constants);
    crdv(m->len_u_veloU_function_constants, m->u_veloU_function_constants);
    crdv(m->len_u_veloL_function_constants, m->u_veloL_function_constants);
    crdv(m->len_u_dcaU_function_constants, m->u_dcaU_function_constants);
    crdv(m->len_u_dcaL_function_constants, m->u_dcaL_function_constants);
    crdv(m->len_lubsource, m->u_lubsource_function_constants);
    /*
     * some more specialized constants (PorousShell)
     */
    crdv(m->len_u_PorousShellClosedPorosity_function_constants,
         m->u_PorousShellClosedPorosity_function_constants);
    crdv(m->len_u_PorousShellClosedRadius_function_constants,
         m->u_PorousShellClosedRadius_function_constants);
    crdv(m->len_u_PorousShellClosedHeight_function_constants,
         m->u_PorousShellClosedHeight_function_constants);
    crdv(m->len_u_PorousShellClosedP0_function_constants,
         m->u_PorousShellClosedP0_function_constants);
    crdv(m->len_u_PorousShellPatm_function_constants, m->u_PorousShellPatm_function_constants);
    crdv(m->len_u_PorousShellPref_function_constants, m->u_PorousShellPref_function_constants);
    crdv(m->len_u_PorousShellCrossKappa_function_constants,
         m->u_PorousShellCrossKappa_function_constants);
    crdv(m->len_u_PorousShellInitPorePres_function_constants,
         m->u_PorousShellInitPorePres_function_constants);
    crdv(m->len_u_PorousShellDiffusivity_function_constants,
         m->u_PorousShellDiffusivity_function_constants);
    crdv(m->len_u_PorousShellRT_function_constants, m->u_PorousShellRT_function_constants);
    crdv(m->len_u_PorousShellHenry_function_constants, m->u_PorousShellHenry_function_constants);
    /*
     * some more specialized constants (thin film)
     */
    crdv(m->len_u_FilmEvap_function_constants, m->u_FilmEvap_function_constants);
    crdv(m->len_u_DisjPress_function_constants, m->u_DisjPress_function_constants);
    crdv(m->len_u_DiffCoeff_function_constants, m->u_DiffCoeff_function_constants);

    crdv(e->len_u_mu, e->u_mu);

    crdv(e->len_u_lambda, e->u_lambda);

    crdv(e->len_u_lame_TempShift, e->u_lame_TempShift);

    crdv(e->len_u_v_mesh_sfs, e->u_v_mesh_sfs);

    crdv(e->len_u_thermal_expansion, e->u_thermal_expansion);

    crdv(e->len_u_solid_viscosity, e->u_solid_viscosity);

    crdv(e->len_u_solid_dil_viscosity, e->u_solid_dil_viscosity);

    e = elc_rs_glob[i];

    crdv(e->len_u_mu, e->u_mu);

    crdv(e->len_u_lambda, e->u_lambda);

    crdv(e->len_u_lame_TempShift, e->u_lame_TempShift);

    crdv(e->len_u_v_mesh_sfs, e->u_v_mesh_sfs);

    crdv(e->len_u_thermal_expansion, e->u_thermal_expansion);

    crdv(e->len_u_solid_viscosity, e->u_solid_viscosity);

    crdv(e->len_u_solid_dil_viscosity, e->u_solid_dil_viscosity);

    crdv(evpl_glob[i]->len_u_plastic_mu, evpl_glob[i]->u_plastic_mu);

    crdv(evpl_glob[i]->len_u_yield, evpl_glob[i]->u_yield);

    v = vn_glob[i];

    crdv(v->len_dg_J_model_wt, v->dg_J_model_wt);

    crdv(v->len_shift, v->shift);

    for (j = 0; j < Num_Var_Init_Mat[i]; j++) {
      crdv(Var_init_mat[i][j].len_u_pars, Var_init_mat[i][j].u_pars);
    }
  }

  // turbulence info
  if (upd->turbulent_info->num_node_sets > 0) {
    ddd_add_member(n, upd->turbulent_info->node_set_ids, upd->turbulent_info->num_node_sets,
                   MPI_INT);
  }

  if (upd->turbulent_info->num_side_sets > 0) {
    ddd_add_member(n, upd->turbulent_info->side_set_ids, upd->turbulent_info->num_side_sets,
                   MPI_INT);
  }

  for (i = 0; i < Num_BC; i++) {
    crdv(BC_Types[i].len_u_BC, BC_Types[i].u_BC);
  }

  for (i = 0; i < num_BC_Tables; i++) {
    crdv(BC_Tables[i]->tablelength, BC_Tables[i]->t);
    crdv(BC_Tables[i]->tablelength, BC_Tables[i]->t2);
    crdv(BC_Tables[i]->tablelength, BC_Tables[i]->f);
  }
  for (i = 0; i < num_MP_Tables; i++) {
    crdv(MP_Tables[i]->tablelength, MP_Tables[i]->t);
    crdv(MP_Tables[i]->tablelength, MP_Tables[i]->t2);
    crdv(MP_Tables[i]->tablelength, MP_Tables[i]->f);
  }
  for (i = 0; i < num_ext_Tables; i++) {
    crdv(ext_Tables[i]->tablelength, ext_Tables[i]->t);
    crdv(ext_Tables[i]->tablelength, ext_Tables[i]->t2);
    crdv(ext_Tables[i]->tablelength, ext_Tables[i]->t3);
    crdv(ext_Tables[i]->tablelength, ext_Tables[i]->f);
  }

  /*   Augmenting condition data tables    */

  for (i = 0; i < nAC; i++) {
    crdv(augc[i].len_AC, augc[i].DataFlt);
    if (augc[i].Aprepro_lib_string_len > 0) {
      ddd_add_member(n, augc[i].Aprepro_lib_string, augc[i].Aprepro_lib_string_len, MPI_CHAR);
    }
  }
  /*   Post-processing parameters    */

  for (i = 0; i < nn_volume; i++) {
    crdv(pp_volume[i]->num_params, pp_volume[i]->params);
  }

  ddd_set_commit(n);

#endif
  return;
}
