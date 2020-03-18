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
 * This routine is the location of stub routines and misc wrappers that
 * Goma needs to communicate with the cpc - chemkin package
 *
 */

/*
 *$Id: mm_chemkin.c,v 5.1 2007-09-18 18:53:42 prschun Exp $
 */

/*************** R O U T I N E S   I N   T H E   F I L E *********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    chemkin_mat_prop_init
 *    ck_decide_vol_chem
 *    chemkin_mat_prop_init
 *    chemkin_not_linked
 *    chemkin_initialize_mp
 ***************************************************************************/

#include <stdio.h>

#include "std.h"
#include "mm_eh.h"

#ifdef USE_CHEMKIN
#include "cpc_defs.h"
#include "ck_chemkin_const.h"
#endif

#include "mm_chemkin.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Globals defined in this file:
 *    
 */

int Chemkin_Needed = FALSE;  /*
			      * If Goma needs the chemkin libraries for
			      * any functionality this flag will be tripped.
			      * It is present whether or not USE_CHEMKIN
			      * ifdef has been activated
			      */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifndef USE_CHEMKIN

void
chemkin_not_linked(char *errString)

    /************************************************************************
     *
     * chemkin_not_linked():
     *
     *    This is an error handler for the case when the use of the
     *    chemkin package is indicated by the input file, but the actual
     *    compile was done without chemkin. A fatal error results
     *
     *    A fatal error exit is triggered
     *
     *    Note: this routine is only compiled into the code when the
     *          USE_CHEMKIN is not defined.
     *
     *  Input
     * -------
     *    errString : string that indicates what  property needed
     *                chemkin
     ************************************************************************/
{
  fprintf(stderr,"ERROR! A chemkin capability has been requested: \n");
  if (errString != NULL) {
    fprintf(stderr, "\t%s\n", errString);
  }
  fprintf(stderr,"\tHowever, CHEMKIN has not been linked in!\n");
  fprintf(stderr,
	  "\tGOMA must be recompiled with the USE_CHEMKIN definition!\n");
  EH(-1, "chemkin not linked in\n");
}
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef USE_CHEMKIN

int
ck_decide_vol_chem(CPC_VOLDOMAIN_STRUCT *volD)

    /*************************************************************************
     *
     * ck_decide_vol_chem():
     *
     *     This routine makes a decision as to what chemkin package will
     * handle calls for the homogeneous chemistry package. Where this routine
     * should reside is debatible. It mixes goma constants and chemkin
     * constants. However, it is really a ck_ routine, because it makes
     * a distinction based on the type of volume domain.
     *
     *  Input
     * ---------
     *  volD   :  Pointer to the current volume domain structure
     *            corresponding to the current material domain.
     *************************************************************************/
{
  CPC_VOLPHASE_STRUCT *volP;
  /*
   *  If this volume domain is anything else than a simple single phase, let
   *  the cpc_ package handle the homogeneous chemistry calls.
   */
  if ((volD->NVolPhase > 1) || (volD->NSurPhase > 0)) return SSM_CHEMKIN_CPC;
  
  /*
   *  Get the single volume phase pointer
   */
  volP = volD->ListVolPhasePtr[0];
  
  /*
   *  If the gas phase boolean is set in the phase structure, then identify
   *  this phase as being a gas phase. Let raw gas phase chemkin handle
   *  the function calls.
   *  (there are many ways to do this!)
   */
  if (volP->GasPhase) return SSM_CHEMKIN_GAS;

  /*
   *  If the liquid boolean is set in the phase structure then identify
   *  this phase as being a liquid phase. Let raw liquid phase chemkin handle
   *  the function calls.
   *  (this is just a placeholder for functionality to be included later)
   */
  if (volP->LiqPhase) return SSM_CHEMKIN_LIQ;

  /*
   *  The default is to use the general package, SSM_CHEMKIN_CPC
   */
  return SSM_CHEMKIN_CPC;
}
#endif 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef USE_CHEMKIN

int
chemkin_mat_prop_init(MATRL_PROP_STRUCT *mat_ptr, int mn,
		      PROBLEM_DESCRIPTION_STRUCT *pd_ptr)
    
    /*************************************************************************
     *
     * chemkin_mat_prop_init:
     *
     *
     *   This routine is called from mm_input_mp.c to initialize a material
     *   structure that has been assigned a chemkin material type. Because
     *   it is called from the input routines, it is only executed on
     *   processor zero.
     *
     *   Currently, this routine initializes the following items in goma's
     *   material structure on proc 0.
     *
     *     MolecularWeightModel[] -> Get the molecular weights from chemkin
     *     molecular_weight[]
     *
     *     SpeciesSourceModel -> Set this variable to indicate that chemkin
     *                           will be called to evaluate source terms.
     *
     *  Input
     * ----------
     *   mat_ptr   : Material Property pointer
     *   mn        : Material number
     *   pd_ptr    : Problem description pointer
     *
     *  Return
     * ----------
     *   1: Successful completion of the routine
     *  -1: error exit
     *************************************************************************/
{
  static int chemkinDatabaseRead = FALSE;
#if 0				/* 991207 pas */
  static int domainMapRead = FALSE;
#endif
  int retn, volD_id, model_id;
  int initGasTransport = TRUE;
  int initLQCHEMKIN = FALSE;
  char *yo = "chemkin_mat_prop_init ERROR";
  CPC_VOLDOMAIN_STRUCT *volD;

  /*
   *  Flip the global flag indicating that chemkin is being used
   */
  Chemkin_Needed = TRUE;

  /*
   *  Read in the cpc data files. These include the chemkin-III data
   *  files and the domain mapping file (currently a virtual entity).
   *  This only needs to be done once, even if there are multiple
   *  domains handled by chemkin. Therefore, we will keep track
   *  of whether we read it by flipping a static boolean,
   *  chemkinDatabaseRead.
   */
  if (!chemkinDatabaseRead) {
    chemkinDatabaseRead = TRUE;
    retn = cpc_initialize_1p(TRUE, ProcID, Num_Proc, initGasTransport,
		 	     initLQCHEMKIN, TRUE, DomainMappingFile);
    if (retn != CPC_SUCCESS) {
      printf("Failure to read chemkin data files: BAIL!\n");
      return -1;
    }
  }

  /*
   *   Next, we have to associate the material ID with the chemkin domain
   *   list. We will do this by taking the character name of the material
   *   stored in the material structure variable, Material_Name, and
   *   then looking up that name in the list of domain names known to
   *   chemkin.
   *      If no match is found, end the code.
   */
  volD_id = cpc_LookupVolD_Name_to_ID(mat_ptr->Material_Name);
  if (volD_id < 0) {
    printf("%s: Mismatch in names!", yo);
    printf("\t\tCouldn't find %s in list of chemkin vol domains\n",
	   mat_ptr->Material_Name);
    EH(-1, yo);
  }

  /*
   *  Now that we have a match for the chemkin domain, we fill up the 
   *  material structure with constituitive models obtained from the
   *  chemkin database.
   */

  /*
   *  Get the volume domain pointer
   */
  volD = cpc_LookupVolD_ID_to_Ptr(volD_id);

  /*
   *  Assign the homogeneous chemistry constituitive model
   */
  model_id = ck_decide_vol_chem(volD);
  cpc_iset(mat_ptr->SpeciesSourceModel, mat_ptr->Num_Species, model_id);

  /*
   *  Copy the species names from the chemkin database for the volume
   *  domain into  the material prop structure.
   */
  mat_ptr->Species_Names = cpc_alloc_VecFixedStrings(mat_ptr->Num_Species,
						     CPC_MAX_NAME_LEN_P1);
  retn = cpc_VD_SpecNames(volD_id, mat_ptr->Species_Names);
  
   /*
    *  Assign the Molecular weights for the species in the volumetric
    *  domain -> Use the default CHEMKIN_MODEL to indicate that the
    *  molecular weights are obtained from chemkin.
    */
  cpc_iset(mat_ptr->MolecularWeightModel, mat_ptr->Num_Species,
	   CHEMKIN_MODEL);
  retn = ck_VD_skwt(volD_id, mat_ptr->molecular_weight);
  if (retn != CPC_SUCCESS) {
    printf("ck_VD_skwt error: %d\n", retn); exit(-1);
  }

  /*
   *  Assign the species variable default and the units for the
   *  species equation. This may change in the future!
   */
  assign_species_var_type( mn, SPECIES_MASS_FRACTION, FALSE);

   /*
    *  Change the defaults on some thermo models from the GOMA
    *  defaults. It's inappropriate to assume that species in
    *  the same slots have the same stoichiometry. Therefore,
    *  all latent heat models and also the vapor pressure model
    *  are inappropriate for CHEMKIN problems. Set these model
    *  types to NO_MODEL.
    */
  cpc_iset(mat_ptr->VaporPressureModel,    mat_ptr->Num_Species,
	   NO_MODEL);
  cpc_iset(mat_ptr->LatentHeatFusionModel, mat_ptr->Num_Species,
	   NO_MODEL);
  cpc_iset(mat_ptr->LatentHeatVapModel,    mat_ptr->Num_Species,
	   NO_MODEL);  
 
  return 1;
}
#endif 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef USE_CHEMKIN

void
chemkin_initialize_mp(void)

    /************************************************************************
     *
     * chemkin_initialize_mp()  
     *
     *       This routine is responsible for copying all of the chemkin
     *    structures from processor 0 to the other processors.
     *    It is currently just a wrapper around the cpc_initialize_mp()
     *    routine. The cpc package can handle this step mostly on its
     *    own.
     *
     *    Error exits from cpc cause fatal program terminations.
     *************************************************************************/
{
  int pl;
  int retn;
  int initGasTransport = TRUE;
  int initLQCHEMKIN = FALSE;
  int initCKExtractGuts = TRUE;
  int infoproc = 0;
#ifdef PARALLEL
#ifdef DEBUG_HKM
  printf("P_%d at barrier before cpc_initialize_mp()\n", ProcID);
  fflush(stdout);
  (void) MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
#ifdef DEBUG_CPC_MP
  pl = 1;
#else
  pl = 0;
#endif
  retn = cpc_initialize_mp(pl, Num_Proc, infoproc, ProcID, initGasTransport,
			   initLQCHEMKIN, initCKExtractGuts);
  if (retn != CPC_SUCCESS) exit(-1);

}
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

