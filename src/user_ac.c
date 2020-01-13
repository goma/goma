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
 

/* $Id: user_ac.c,v 5.4 2010-07-21 16:39:27 hkmoffa Exp $ */

/* GOMA include files */


#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_bc_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mpi.h"
#include "rf_bc.h"
#include "user_ac.h"

#define GOMA_USER_AC_C

/****************************************************************************/
/*

  ADAPTED FROM user_bc.c

  BY IAN GATES
  2/98 - 9/98

*/

void user_aug_cond_residuals(int iAC, 
			     double *x, 
                             double *xdot,
                             double delta_t,
                             double time_value,
                             double **x_sens_p,
			     double *AC, 
			     int *have_bAC, int *have_cAC, int *have_dAC,
			     double **bAC, double **cAC, double **dAC,
			     Exo_DB *exo,
			     Dpi *dpi,
			     Comm_Ex *cx)
{
  /*  dbl inventory, target; */

  
  /* Goma Users and Developers: 


      Add your ACs in the section below the double $$$$$$$ line.


     ==========================
     Miscellaneous information: 
        How to extract values of solution variables from database
	  ndof = 0;
	  i = Index_Solution(100, MESH_DISPLACEMENT1, 0, ndof, -1, matrix ID);
	          or for a specific material id:
          i = Index_Solution(100, MESH_DISPLACEMENT1, 0, ndof, mat_id, matrix ID);  
     ==========================

      Please read the following 2 notes and 1 example first!!                   */

  /******************************************************************************/

  /* Note 1 */
  /* The primary control during application of ACs is through the TYPE of AC
     you have specified in the Goma input file, i.e., AC = BC/MT/VC/FC.

     The AC.Type is read from the input and calculation flow is directed
     accordingly through the if-else-endif structure below.  The user
     inserts his/her ACs between the if block curly braces { } in the
     appropriate branch of the if-else-endif structure below according
     to whether the AC is a BC, MT, VC or FC. If multiple types are present,
     each type is inserted in its appropriate branch; if multiple ACs of
     the same type are defined, they are inserted in the same branch. Note, in
     this case, the float_list can be used to set a sub-model as shown by the
     example immediately following this comment.)                               */

  /******************************************************************************/
  /* An example of multiple BC_Type ACs:                                        */
  /******************************************************************************/
  /*                                                                            */
  /**  loop over all AC conditions **/
  /* for(iAC=0;iAC<nAC;iAC++)                                                   */
  /*   {                                                                        */
  /*                                                                            */
  /* BC augmenting condition */
  /*     if (augc[iAC].Type == AC_USERBC )                                              */
  /*     {                                                                      */
  /*                                                                            */
  /*      Set submodel id (first value in AC float_list)                        */
  /*        model_id = 0;                                                       */
  /*        if(augc[iAC].len_AC > 0 ) model_id = (int)augc[iAC].DataFlt[0];     */
  /*                                                                            */
  /*      Submodel 1                                                            */
  /*        if(model_id == 1)                                                   */
  /*          {                                                                 */
  /*           int nsp,ns_id,bc_id,cabc_id;                                     */
  /*           dbl radius, surf_tens;                                           */
  /*                                                                            */
  /*           ns_id = (int) augc[iAC].DataFlt[1];                              */
  /*           bc_id = (int) augc[iAC].DataFlt[2];                              */
  /*           cabc_id = (int) augc[iAC].DataFlt[3];                            */
  /*           surf_tens = BC_Types[cabc_id].BC_Data_Float[0];                  */
  /*           nsp       = match_nsid(ns_id);                                   */
  /*           k         = Proc_NS_List[Proc_NS_Pointers[nsp]];                 */
  /*                                                                            */
  /*           for (j = 0; j < Proc_NS_Count[nsp]; j++)                         */
  /*             {                                                              */
  /*               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];                   */
  /*               i = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);        */
  /*               EH(i, "Could not resolve Index_Solution.");                  */
  /*               radius             = Coor[1][k] + x[i];                      */
  /*             }                                                              */
  /*                                                                            */
  /*           AC[iAC] = BC_Types[bc_id].BC_Data_Float[0]-surf_tens/radius;     */
  /*          }                                                                 */
  /*                                                                            */
  /*      NB:the following AC can now be handled using the FC type AC.          */
  /*         It is included here as an example of changing a BC in one location */
  /*         to control a flux (force) on another boundary.                     */
  /*      Submodel 2 (flux condition)                                           */
  /*        if(model_id == 2)                                                   */
  /*          {                                                                 */
  /*           int ss_id,flux_type, blk_id, species_num;                        */
  /*           double n_load,n_load1;                                           */
  /*           char *filenm;                                                    */
  /*                                                                            */
  /*           ss_id=(int) augc[iAC].DataFlt[1];                                */
  /*           flux_type= (int) augc[iAC].DataFlt[2];                           */
  /*           blk_id= (int) augc[iAC].DataFlt[3];                              */
  /*           n_load1 = augc[iAC].DataFlt[4];                                  */
  /*           species_num=0;                                                   */
  /*                                                                            */
  /*           for(i=0;i<NumUnknowns[pg->imtrx];i++){                           */
  /*             cAC[0][i]=0.0 ;         }                                      */
  /*                                                                            */
  /*           af->Assemble_Jacobian = TRUE;                                    */
  /*           n_load = evaluate_flux(exo, dpi,  ss_id, flux_type,              */
  /*                                  NULL, blk_id, species_nu                  */
  /*                        ,filenm,x,xdot,&cAC[0][0], delta_t,time_value,0);   */
  /*           AC[iAC] = n_load - n_load1;                                      */
  /*                                                                            */
  /*           *have_cAC = TRUE;                                                */
  /*          }                                                                 */
  /* MT augmenting condition                                                    */
  /*   else if (augc[iAC].Type == AC_USERMAT )                                  */
  /*     {                                                                      */
  /*     }                                                                      */
  /* VC augmenting condition                                                    */
  /*   else if (augc[iAC].Type == AC_VOLUME )                                   */
  /*     {                                                                      */
  /*       inventory = augc[iAC].evol;                                          */
  /*       target = augc[iAC].CONST;                                            */
  /*       AC[iAC] =  -target + inventory;                                      */
  /*     }                                                                      */
  /*                                                                            */
  /*   }  End of Loop on ACs                                                    */
  /******************************************************************************/
  /* End of Example */ 
  /******************************************************************************/


  /* Note 2 */
  /* The BC Augmenting Condition implemented below is used by three (3) of the
     Goma Test Suite Problems. Please reinstate this AC before you check-in a
     user_ac.c.                                                                 */

  /* Please DO NOT CHECK IN A CHANGED user_ac.c  */
  /* unless you have added some necessary infrastructure changes. */

  /******************************************************************************/
  /* End of Notes; Source Code follows. */ 
  /******************************************************************************/


  /* User ACs should be added in the section below the next 2 lines of $$$$$$   */

  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

  /* Comment-out next line if 1 or more ACs is implemented; uncomment if no ACs  */
  /*EH(-1,"You have provided no augmenting condition in user_ac.c"); */

     /* NB, the volume constraint residual is based on an integrated quantity
        of either volume, overall mass, or component mass; hence, it is not
        associated with the x[] vector explicitly and has to be constructed
        via the AC_Information struct. ACSun Feb, 1999. Example:             */


  /**  loop over all AC conditions **/

  if (augc[iAC].Type == AC_USERBC )           /* BC augmenting condition */
    {
      AC[iAC] = BC_Types[4].BC_Data_Float[0]-2.0;  /* Goma Test Suite */
    }
  else if (augc[iAC].Type == AC_USERMAT )      /* MT augmenting condition */
    {
    }

 /* VC augmenting condition */

  /*  else if (augc[iAC].Type == AC_VOLUME )     
    {
      inventory = augc[iAC].evol;
      target = augc[iAC].CONST; 
      AC[iAC] =  -target + inventory; 
      } */
  /*
  else if (augc[iAC].Type == AC_FLUX )      
    {
      int mn = map_mat_id( augc[iAC].MTID );

      nAC = 1;

      inventory = evaluate_flux( exo,
                                 dpi, 
				 augc[iAC].SSID,
				 augc[iAC].MFID,
				 NULL,
				 mn,
				 augc[iAC].COMPID,
				 NULL,
				 x,
				 xdot,
				 (*have_cAC == TRUE ? NULL : cAC[iAC]),
				 delta_t,
				 time_value,
				 0);
	  
      target = augc[iAC].CONST; 
      AC[iAC] =  -target + inventory; 

      if( (iAC + 1) == nAC )     *have_cAC = TRUE;

    }
  */

} /* END of routine user_aug_cond_residuals                                 */
/****************************************************************************/

/*
 *  user_aug_cond_volume_residuals():
 *
 *    This routine implements an expansion of the Volume Constraint
 *    augmented condition.  Using this routine, you can implement any
 *    function that employs the total volume of a material calculated 
 *    using the VC condition.
 *
 *    To access this implementation you need to specify a Volume
 *    Constraint type of 11, 12, or 13, instead of the types 1, 2, or 3.
 * 
 *    see pg. 25 of Advanced" Capabilities in Goma 5.0 - Augmenting
 *            Conditions, automatic Continuation and Linear Stability Analysis,"
 *            SAND 2--6-7304
 *
 *    The total integrated volume is fed into the routine via the
 *    parameter  augc[iAC].evol. The derivative of the total volume
 *    with respect to the soln unknowns are located in augc[iAC].d_evol_dx[i].
 *   augc[iAC].CONSTV contains the constant supplied on the AC = VC card.
 *
 *  parameters
 *
 * Input:
 *        iAC      index of the augmented condition
 *        x[]      Raw solution vector
 *        xdot[]   Raw solution time derivative vector
 *        deltat   delta t
 *        time_value  Time 
 *        numProcUnknowns  Number of unknowns in the soln vector
 *                       on this processor
 *        exo            Exodus file
 *        dpi            helper struct
 *        cx             helper struct
 *       
 * Output:
 *        AC[iAC]        Value of the residual
 *        cAC[iAC][i]    Jacobian entries for the dependence of aug residual, iAC
 *                       on the solution unknown, i
 *        dAC[iAC][jAC]  Jacobian entries for the dependence of aug residual, iAC
 *                       on the aug unknown, jAC 
 *
 *  The default implementation in this routine just duplicates the
 *  volume type 1, 2, and 3 VC conditions.
 */
void user_aug_cond_volume_residuals(const int iAC, 
				    const double * const x, 
				    const double * const xdot,
				    const double delta_t,
				    const double time_value,
				    const double * const x_AC, 
				    double * const AC,
				    double ** const cAC,
				    double **const dAC,
				    const int numProcUnknowns,
				    const Exo_DB * const exo,
				    const Dpi * const dpi,
				    const Comm_Ex * const cx)
{
  /*********************************************************************/
  int i, jAC;
  double inventory, target;
#ifdef PARALLEL
  double global_inventory;
#endif
  inventory = augc[iAC].evol;
#ifdef PARALLEL
  if (Num_Proc > 1) 
    {
      MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
      inventory = global_inventory;
    }
#endif
  for (i = 0; i < numProcUnknowns; i++)
    {
      cAC[iAC][i] = augc[iAC].d_evol_dx[i];
    }
  
  for (jAC = 0; jAC < nAC; jAC++)
    {
      dAC[iAC][jAC] = 0.0;
    }
      
  target = augc[iAC].CONSTV; 
  AC[iAC] =  -target + inventory; 
}
/****************************************************************************/
