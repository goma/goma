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

/* $Id: user_ac.c,v 5.4 2010-07-21 16:39:27 hkmoffa Exp $ */

/* GOMA include files */

#include "user_ac.h"
#include "dp_types.h"
#include "dpi.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mpi.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#include "std.h"
#include "mm_augc_util.h"
#include "el_geom.h"
#include "mm_mp.h"
#include "wr_side_data.h"

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
  dbl target;
  int i, j, k;
  int model_id;


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
  /*               GOMA_EH(i, "Could not resolve Index_Solution.");                  */
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
  /*GOMA_EH(GOMA_ERROR,"You have provided no augmenting condition in user_ac.c"); */

  /* NB, the volume constraint residual is based on an integrated quantity
     of either volume, overall mass, or component mass; hence, it is not
     associated with the x[] vector explicitly and has to be constructed
     via the AC_Information struct. ACSun Feb, 1999. Example:             */

  /**  loop over all AC conditions **/

  if (augc[iAC].Type == AC_USERBC )           /* BC augmenting condition */
    {
        model_id = 0;
        if(augc[iAC].len_AC > 0 )
		{
		if( augc[iAC].BCID == APREPRO_AC_BCID)
			{
 			model_id = (int)augc[iAC].DataFlt[1];
			}
			else	
			{
 			model_id = (int)augc[iAC].DataFlt[0];
			}	
		}
	if(model_id == 0)
	{
      AC[iAC] = BC_Types[4].BC_Data_Float[0]-2.0;  /* Goma Test Suite */
	}

	else if (model_id == 1)    /* This is for ePQR problem to find the location of ns_id x-location 
                                      that gives a given coating_thicknes (coat_th=ytop-ybottom).
                                      x2=x1 along a line segmant and  y2=y1+coat_th at the end point 
                                      ns_id2 follows ns_id1
                                      ytop=yt lies on the circle (x-Xc)^2+(y-Yc)^2=R^2
                                      yt=Yc-sqrt(R^2-(xt-Xc)^2) 
                                      yb=ybottom : y at the bottom roll (rubber roll)
                                      yt-yb=coat_th=Yc-yb -sqrt(R^2-(xt-Xc)^2)
                                      which gives xt as (take + sign to be on the right of center line between rolls)
                                      x=Xc+sqrt(R^2-(yt-Yc)^2)
                                      x=Xc+sqrt(R^2-(yb+coat_th-Yc)^2)
                                      dx=x-X=Xc-X+sqrt(R^2-(yb+coat_th-Yc)^2)
                                      */
	{
          int nsp,ns_id1,ns_id2;
          double coat_th, R, Xc, Yc, X, xb, yb, xt;
          ns_id1 = (int) augc[iAC].DataFlt[1];
          ns_id2 = (int) augc[iAC].DataFlt[2];
          coat_th = augc[iAC].DataFlt[3];
          R = augc[iAC].DataFlt[4];
          Xc = augc[iAC].DataFlt[5];
          Yc = augc[iAC].DataFlt[6];

          nsp = match_nsid(ns_id1);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               X = Coor[0][k];
               xb = Coor[0][k] + x[i];
               yb = Coor[1][k] + x[i+1];
             }

          nsp = match_nsid(ns_id2);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xt = Coor[0][k] + x[i];
             }

          if (R*R > SQUARE(yb+coat_th-Yc))
            {
              AC[iAC] = xt - (Xc+sqrt(R*R-SQUARE(yb+coat_th-Yc)));
            }
          else
            {
              AC[iAC] = xt - Xc;
            }
	}

        else if(model_id == 2)                 
        {                               
          int ss_id,flux_type, blk_id1, blk_id2, species_num; 
          double flux, flux1, flux2;                   
                                                  
          ss_id=(int) augc[iAC].DataFlt[1];     
          flux_type= (int) augc[iAC].DataFlt[2];
          blk_id1= (int) augc[iAC].DataFlt[3];  
          blk_id2= (int) augc[iAC].DataFlt[4];  
          flux = augc[iAC].DataFlt[5];     
          species_num=0;                     
                                             
          for(i=0;i<NumUnknowns[pg->imtrx];i++){      
          cAC[0][i]=0.0 ;         }     
                                         
          af->Assemble_Jacobian = TRUE; 
          flux1 = evaluate_flux(exo, dpi,  ss_id, flux_type,
                          NULL, blk_id1, species_num  
                    ,NULL,FALSE,x,xdot,&cAC[0][0], delta_t,time_value,0);
          flux2 = evaluate_flux(exo, dpi,  ss_id, flux_type,
                          NULL, blk_id2, species_num  
                    ,NULL,FALSE,x,xdot,&cAC[0][0], delta_t,time_value,0);
          AC[iAC] = flux1 + flux2 - flux; 
          *have_cAC = TRUE;         
        }                         

	else if (model_id == 3)    /* ns_id2 follows ns_id1 in the x-direction (x2=x1+dx) */
	{
          int nsp,ns_id;
          double xpoint,dx;
          ns_id = (int) augc[iAC].DataFlt[1];

          nsp = match_nsid(ns_id);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               dx = x[i];
               xpoint = Coor[0][k] + x[i];
             }

          AC[iAC] = dx;
	}

	else if (model_id == 4)    /* ns_id2 follows ns_id1 in the idir-direction */
	{
 	  int nsp, ns_id1, ns_id2, idir, index, pt_count1=0, pt_count2=0;
          double coord1_sum=0, coord2_sum=0, dy=0, penalty=1./BIG_PENALTY;;
          ns_id1 = (int) augc[iAC].DataFlt[1];
          ns_id2 = (int) augc[iAC].DataFlt[2];
          idir = (int) augc[iAC].DataFlt[3];
	  if(augc[iAC].len_AC > 3) dy = augc[iAC].DataFlt[4];
	  if(augc[iAC].len_AC > 4) penalty = augc[iAC].DataFlt[5];

          nsp = match_nsid(ns_id1);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1+idir-1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
	       coord1_sum += Coor[idir-1][k] + x[i];
	       pt_count1++;
             }

          nsp = match_nsid(ns_id2);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1+idir-1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
	       coord2_sum += Coor[idir-1][k] + x[i];
	       pt_count2++;
             }
          *have_cAC = TRUE;
          for(i=0;i<NumUnknowns[pg->imtrx];i++){    cAC[iAC][i] = 0.0; } 
          nsp = match_nsid(ns_id2);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1+idir-1, 0, 0, -1, pg->imtrx);
               cAC[iAC][i] = penalty*pt_count1;
             }

          AC[iAC] = pt_count1*coord2_sum - pt_count2*(coord1_sum + dy*pt_count1);
          AC[iAC] *= penalty;
/*fprintf(stderr,"AC %d %g %g %g %g %d %d\n",iAC,coord1_sum/pt_count1,coord2_sum/pt_count2,
          coord2_sum/pt_count2-coord1_sum/pt_count1,AC[iAC],pt_count1,pt_count2);*/
	}
  
	else if (model_id == 5)    /* ns_id2 follow ns_id = ns_id1: old implementation*/
	{
          int nsp,ns_id,ns_id2;
          double xpoint,xpoint1;
          ns_id = (int) augc[iAC].DataFlt[1];
          ns_id2 = (int) augc[iAC].DataFlt[2];
          nsp            = match_nsid(ns_id);
          k           = Proc_NS_List[Proc_NS_Pointers[nsp]];

          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint = Coor[0][k] + x[i];
             }

          nsp = match_nsid(ns_id2);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];

          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint1 = Coor[0][k] + x[i];
             }

          AC[iAC] = xpoint1-xpoint;
	}
  
        else if (model_id == 6)
                {
                dbl temp, spec_flux, heat_flux;
                double moles_init, mass_cum, delta_mass, theta;
		double termK, rad_B, Pamb=1.0, Rgas=82.07, constK;
		double cgs=1013250.0,surf_tens, Tref=298.15;
                int ss_id, blk_id, ns_id, nsp, T_index, extra_bc;
                rad_B = augc[iAC].DataFlt[0];
                ss_id = (int) augc[iAC].DataFlt[3];
                blk_id = (int) augc[iAC].DataFlt[4];
                moles_init = augc[iAC].DataFlt[2];
                mass_cum = augc[iAC].DataFlt[5];
                theta = augc[iAC].DataFlt[7];
		extra_bc = augc[iAC].DFID;
		constK = 6./(1+cos(theta))/(2.-cos(theta));
                temp = augc[iAC+1].DataFlt[0];
		surf_tens = mp->u_surface_tension[0]+
			mp->u_surface_tension[1]*(temp-Tref);
                for(i=0;i<NumUnknowns[pg->imtrx];i++)     {cAC[iAC][i]=0.0 ;}   
                for(i=0;i<NumUnknowns[pg->imtrx];i++)     {bAC[iAC][i]=0.0 ;}   
                af->Assemble_Jacobian = TRUE;  
                spec_flux = evaluate_flux(exo, dpi, ss_id, 7, "SPECIES_FLUX", 
                                blk_id, 0 ,NULL,0,x,xdot,&cAC[iAC][0], delta_t,
                                time_value,0); 
                heat_flux = evaluate_flux(exo, dpi, ss_id, 8, "HEAT_FLUX", 
                                blk_id, 0 ,NULL,0,x,xdot,&bAC[iAC][0], delta_t,
                                time_value,0); 
		delta_mass = spec_flux*delta_t;
		termK = Rgas*(moles_init+(mass_cum+delta_mass)
                              /mp->molecular_weight[0]);

                AC[iAC] = 2.*M_PIE*(1.+cos(theta))*
                         (Pamb*CUBE(rad_B)+ surf_tens/cgs*SQUARE(rad_B))
                          /constK-termK*temp;
                AC[iAC] *= 1.0E+09;
  
                dAC[iAC][iAC]=2.*M_PIE*(1.+cos(theta))*(3.*SQUARE(rad_B)*Pamb+
			surf_tens/cgs*2.*rad_B)/constK*1.0E+09;
                dAC[iAC][iAC+1]=-termK*1.0E+09;
                *have_dAC = TRUE;
  
                for(i=0;i<NumUnknowns[pg->imtrx];i++)
                      {cAC[iAC][i] *= 
			(-Rgas/mp->molecular_weight[0]*temp*delta_t)*1.0E+09;}   
                *have_cAC = TRUE;   
fprintf(stderr,"AC1 %g %g %g %g %g\n",rad_B,delta_mass,termK,temp,AC[iAC]);

                for(i=0;i<NumUnknowns[pg->imtrx];i++){    bAC[iAC][i] = 0.0; }   
                *have_bAC = FALSE;   
                augc[iAC].DataFlt[6] = delta_mass;
                }
        else if (model_id == 61)
                {
                dbl temp_init, temp, spec_flux,moles_init;
                double mass_cum, delta_heat, heat_cum, heat_cap;
		double Pamb=1.0, Rgas=82.07;
		double alpha_air=5./2., alpha_sol=7./2.;
		double cgs=1013250.0, heat_flux;
                int ss_id, blk_id, extra_bc;
                temp = augc[iAC].DataFlt[0];
                temp_init = augc[iAC].DataFlt[2];
                ss_id = (int) augc[iAC].DataFlt[3];
                blk_id = (int) augc[iAC].DataFlt[4];
                heat_cum = augc[iAC].DataFlt[5];
                mass_cum = augc[iAC-1].DataFlt[5];
                heat_cap = augc[iAC].DataFlt[7];
                moles_init = augc[iAC].DataFlt[8];
		extra_bc = augc[iAC].DFID;
                for(i=0;i<NumUnknowns[pg->imtrx];i++)     {cAC[iAC][i]=0.0 ;}   
                for(i=0;i<NumUnknowns[pg->imtrx];i++)     {bAC[iAC][i]=0.0 ;}   
                af->Assemble_Jacobian = TRUE;  
                spec_flux = evaluate_flux(exo, dpi, ss_id, 7, "SPECIES_FLUX", 
                                blk_id, 0 ,NULL,0,x,xdot,&bAC[iAC][0], delta_t,
                                time_value,0); 
                heat_flux = evaluate_flux(exo, dpi, ss_id, 8, "HEAT_FLUX", 
                                blk_id, 0 ,NULL,0,x,xdot,&cAC[iAC][0], delta_t,
                                time_value,0); 
  
		delta_heat = (spec_flux*mp->latent_heat_vap[0]+heat_flux)*delta_t;

                AC[iAC] = heat_cum+delta_heat
                           -(alpha_air*moles_init+alpha_sol*
                        (mass_cum+spec_flux*delta_t)/mp->molecular_weight[0])
                         *Rgas*(temp-temp_init);
                AC[iAC] *= 1.0E+09;
  
                dAC[iAC][iAC]=-(alpha_air*moles_init+alpha_sol*
                   (mass_cum+spec_flux*delta_t)/mp->molecular_weight[0])*Rgas;
                dAC[iAC][iAC] *= 1.0E+09;
fprintf(stderr,"AC2 %g %g %g %g %g\n",temp,temp_init,delta_heat,spec_flux,AC[iAC]);
                *have_dAC = TRUE;
  
                for(i=0;i<NumUnknowns[pg->imtrx];i++)
                      {cAC[iAC][i] *= 1.0E+09*delta_t;}   
                for(i=0;i<NumUnknowns[pg->imtrx];i++)
                      {cAC[iAC][i] += 1.0E+09*mp->latent_heat_vap[0]*bAC[iAC][i]*delta_t;}   
                for(i=0;i<NumUnknowns[pg->imtrx];i++)
                      {cAC[iAC][i] -= 1.0E+09*alpha_sol*delta_t/mp->molecular_weight[0]
				*bAC[iAC][i]*Rgas*(temp-temp_init);}   
                *have_cAC = TRUE;   

                for(i=0;i<NumUnknowns[pg->imtrx];i++){    bAC[iAC][i] = 0.0; }   
                *have_bAC = FALSE;   
                augc[iAC].DataFlt[6] = delta_heat;
                }
	else if (model_id == 7)
 	{
 	int nsp, ns_id, idir, index, pt_count;
        double xpoint, coord, coord_mean;
        ns_id = (int) augc[iAC].DataFlt[1];
        idir = (int) augc[iAC].DataFlt[2];
        coord = augc[iAC].DataFlt[3];
      	nsp = match_nsid(ns_id);
      	k = Proc_NS_List[Proc_NS_Pointers[nsp]];
	pt_count = 0;
	coord_mean = 0.0;
        for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
	       index = MESH_DISPLACEMENT1-1+idir;
               i = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint = Coor[idir-1][k] + x[i];
	       coord_mean += xpoint;
	       pt_count++;
/*fprintf(stderr,"AC %d %d %g\n",j,i,xpoint); */
              }
  	AC[iAC] = coord_mean - pt_count*coord;
/*fprintf(stderr,"AC %g %d %g %g\n",coord_mean,pt_count,coord,AC[iAC]);*/
 	}
 	else if (model_id == 8)
 	{
 	int nsp,ns_id,ns_id2,index;
        double sca,dca,die_a;
 	double yscl,ydcl,rad;
        ns_id = (int) augc[iAC].DataFlt[2];
        ns_id2 = (int) augc[iAC].DataFlt[3];
        rad = augc[iAC].DataFlt[4];
        sca = augc[iAC].DataFlt[5];
        dca = augc[iAC].DataFlt[6];
        die_a = augc[iAC].DataFlt[7];
 
      	nsp = match_nsid(ns_id);
      	k = Proc_NS_List[Proc_NS_Pointers[nsp]];
 
         for (j = 0; j < Proc_NS_Count[nsp]; j++)
              {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
 	       index = MESH_DISPLACEMENT2;
               i = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               yscl = Coor[1][k] + x[i];
              }
     	nsp = match_nsid(ns_id2);
      	k = Proc_NS_List[Proc_NS_Pointers[nsp]];
 
        for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
 	       index = MESH_DISPLACEMENT2;
               i = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               ydcl = Coor[1][k] + x[i];
              }
   	AC[iAC] = yscl-ydcl+rad*(cos(sca-die_a)+cos(dca));
	}
	else if (model_id == 9)
  	{
  	int nsp,ns_id,idx,idy,N;
 	int id[13]={12,11,10,9,8,7,6,5,4,3,2,1,0};
 	double xpt[13],ypt[13],cee,phi[3];
        double xpoint,ypoint,delta,point;
 	double roll_rad,inter;
        ns_id = (int) augc[iAC].DataFlt[1];
        roll_rad =  augc[iAC].DataFlt[2];
        inter =  augc[iAC].DataFlt[3];
        point =  augc[iAC].DataFlt[4];
       	nsp = match_nsid(ns_id);
       	k = Proc_NS_List[Proc_NS_Pointers[nsp]];
 	N = Proc_NS_Count[nsp];
         for (j = 0; j < N; j++)
              {
                k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
                idx = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
                idy = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
                GOMA_EH(idx, "Could not resolve index_solution.");
                GOMA_EH(idy, "Could not resolve index_solution.");
                xpoint = Coor[0][k] + x[idx];
                ypoint = Coor[1][k] + x[idy];
 		xpt[id[j]] = xpoint;
 		ypt[id[j]] = ypoint;
 		delta= ypoint-(roll_rad-inter-sqrt(roll_rad*roll_rad-point*point));
               }
 	for( i=0; i < N ; i+=2)
         {
           cee = (point-xpt[i])/(xpt[i+2] - xpt[i]);
           if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
               || (cee > 1.0 && i == N-2))
           {
             phi[0]=2.*cee*cee-3.*cee+1.;
             phi[1]=-4.*cee*cee+4.*cee;
             phi[2]=2.*cee*cee-cee; 
             ypoint = ypt[i]*phi[0]+ypt[i+1]*phi[1]+ypt[i+2]*phi[2];
             break;
           }
         }
 
 		delta= ypoint-(roll_rad-inter-sqrt(roll_rad*roll_rad-point*point));
 
   	AC[iAC] = delta;
  	}
        /* Compression of trapped air-  RSBay Feb14-06 based on input from SJendoubi  */
        /*  p V^gamma = p0 V0^gamma   */
        else if (model_id == 10)
          {
           dbl Press0, Volume0, VolumeAdd, gamma, volume, press;
           int blk_id;
           Press0 = augc[iAC].DataFlt[1];
           Volume0  = augc[iAC].DataFlt[2];
           VolumeAdd  = augc[iAC].DataFlt[3];
           blk_id = (int) augc[iAC].DataFlt[4];
           gamma  = augc[iAC].DataFlt[5];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];
  
              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;
           volume = evaluate_volume_integral(exo, dpi, 8, "POSITIVE_FILL", blk_id, 0
                             ,NULL,NULL,0,&cAC[iAC][0],x,xdot, delta_t,time_value,0)-VolumeAdd;
  
           AC[iAC] = press*pow(volume,gamma) - Press0*pow(Volume0,gamma);

           dAC[iAC][iAC]=pow(volume,gamma);
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
                   *have_bAC = TRUE;
          }
        /* Simple Nip Pressure Profile-  RSBay Feb15-06                    */
        /*      Used to augment pressure boundary conditions               */
        /*  pressure increases linearly from P0 to Pmax in time t1         */
        /*  then decreases back again linearly to time t2 and holds at P0  */
        else if (model_id == 11)
          {
           dbl Press0, PressMax, t1, t2, press;
           Press0 = augc[iAC].DataFlt[1];
           PressMax  = augc[iAC].DataFlt[2];
           t1 = augc[iAC].DataFlt[3];
           t2 = augc[iAC].DataFlt[4];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;

                       if (time_value < t1)
             {
             AC[iAC] = press -Press0 -(PressMax-Press0)*time_value/t1;
             }
                       else if (time_value < t2)
                         {
             AC[iAC] = press -PressMax +(PressMax-Press0)*(time_value-t1)/(t2-t1);
                         }
                       else
                         {
             AC[iAC] = press - Press0;
                         }
           dAC[iAC][iAC]=1.0;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
             *have_bAC = TRUE;
          }
        /* Simple Nip Force Profile-  RSBay Mar13-06                       */
        /*     Used to augment velocity boundary conditions                */
        /*  Force over area increases linearly from P0 to Pmax in time t1  */
        /*  then decreases back again linearly to time t2 and holds at P0  */
        else if (model_id == 12)
          {
           dbl Press0, PressMax, t1, t2, area, force;
           int blk_id, ss_id;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           Press0 = augc[iAC].DataFlt[3];
           PressMax  = augc[iAC].DataFlt[4];
           t1 = augc[iAC].DataFlt[5];
           t2 = augc[iAC].DataFlt[6];
           area = augc[iAC].DataFlt[7];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

           force = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

             *have_cAC = TRUE;

                       if (time_value < t1)
                         {
             target = area*( Press0 +(PressMax-Press0)*time_value/t1);
                         }
                       else if (time_value < t2)
                         {
             target = area*(PressMax -(PressMax-Press0)*(time_value-t1)/(t2-t1));
                         }
                       else
                         {
             target = area*Press0;
			 }
             AC[iAC] =  -target + force; 
	  }
        /* Arbitrary Force Profile-  RSBay Feb2-11                            */
        /*     Used to augment velocity boundary conditions                   */
        /*  Force over area follows table input with interpolation            */
	/*  Below t1 it holds at P1, and above t2 it holds at P2              */
        else if (model_id == 13)
          {
           dbl Press1, Press2, t1, t2, area, force;
           int blk_id, ss_id;
              /* New line for table input */
              double slope, table_abs[1], pressure_target;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           Press1 = augc[iAC].DataFlt[3];
           Press2  = augc[iAC].DataFlt[4];
           t1 = augc[iAC].DataFlt[5];
           t2 = augc[iAC].DataFlt[6];
           area = augc[iAC].DataFlt[7];
              /* New lines for table input */
	      table_abs[0]=time_value;
	      pressure_target = interpolate_table( AC_Tables[0], table_abs, &slope, NULL);
  
	   for(i=0;i<NumUnknowns[pg->imtrx];i++){      
	       cAC[iAC][i]=0.0 ; }   
               af->Assemble_Jacobian = TRUE; 
                                     
           force = evaluate_flux( exo,
                                  dpi, 
		  		  ss_id,
				  3,
				  "FORCE_X",
				  blk_id,
				  0,
				  NULL,
				  FALSE,
				  x,
				  xdot,
				  &cAC[iAC][0],
				  delta_t,
				  time_value,
				  0);

	     *have_cAC = TRUE; 

		       if (time_value < t1)
			 {
             target = area*Press1;
			 }
		       else if (time_value < t2)
			 {
             target = area*pressure_target;
			 }
		       else
			 {
             target = area*Press2;
			 }
             AC[iAC] =  -target + force; 
	  }
/*-------------------------------------------------------------------------*/
        /* computes the pressure based on p=P0*V0/V where                  */
        /* Volume= Volume0-(xpoint-x0)*area;                               */
        else if(model_id == 14)                 /* Slah user_ac */
        {
          double x0, area, Press0, Press, Volume0, Volume;
          int nsp,ns_id,idir,index;
          double xpoint;
          x0 = augc[iAC].DataFlt[1];
          area  = augc[iAC].DataFlt[2];
          Press0 = augc[iAC].DataFlt[3];
          Volume0 = augc[iAC].DataFlt[4];
          ns_id = (int) augc[iAC].DataFlt[5];
          idir = (int) augc[iAC].DataFlt[6];


           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

          nsp = match_nsid(ns_id);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
               {
                 k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
                 index = MESH_DISPLACEMENT1-1+idir;
                 i = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
                 GOMA_EH(i, "Could not resolve index_solution.");
                 xpoint = Coor[idir-1][k] + x[i];
                }

          Volume= Volume0-(xpoint-x0)*area;

          Press =Volume0*Press0/Volume;
          AC[iAC] = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID] - Press;
        }
/*-------------------------------------------------------------------------*/
        /* Compression of trapped air-     slah */
        /*  p V^gamma = p0 V0^gamma             */
        /* V = v_air = volume1 - v_liquid       */
        else if (model_id == 15)
          {
           dbl Press0, Volume0, Volume1, gamma, volume, press, v_liquid;
           int blk_id;
           Press0 = augc[iAC].DataFlt[1];
           Volume0  = augc[iAC].DataFlt[2];
           Volume1  = augc[iAC].DataFlt[3];
           blk_id = (int) augc[iAC].DataFlt[4];
           gamma  = augc[iAC].DataFlt[5];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;
           /*volume = evaluate_volume_integral(exo, dpi, 8, "POSITIVE_FILL", blk_id, 0*/
           v_liquid = evaluate_volume_integral(exo, dpi, 0, "VOLUME", blk_id, 0
                             ,NULL,NULL,0,&cAC[iAC][0],x,xdot, delta_t,time_value,0);

           volume = Volume1 - v_liquid;

           AC[iAC] = press*pow(volume,gamma) - Press0*pow(Volume0,gamma);

           dAC[iAC][iAC]=pow(volume,gamma);
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
                   *have_bAC = TRUE;
          }

/*-------------------------------------------------------------------------*/
        /* Pressure profile                 slah jendoubi                  */
        /*      Used to augment pressure boundary conditions               */
        /*  pressure increases linearly from P0 to Pmax in time t1         */
        /*  then holds at Pmax                                              */
        else if (model_id == 16)
          {
           dbl Press0, PressMax, t1, t2, press;
           Press0 = augc[iAC].DataFlt[1];
           PressMax  = augc[iAC].DataFlt[2];
           t1 = augc[iAC].DataFlt[3];
           t2 = augc[iAC].DataFlt[4];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;

                       if (time_value < t1)
             {
             AC[iAC] = press -Press0 -(PressMax-Press0)*time_value/t1;
             }
                       else
                         {
             AC[iAC] = press - PressMax;
                         }
           dAC[iAC][iAC]=1.0;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
             *have_bAC = TRUE;
          }
/*-------------------------------------------------------------------------*/
        /* Pressure profile                 slah jendoubi                  */
        /*      Used to augment pressure boundary conditions               */
        /*  pressure is given by Press0*(1-exp(-time/t1)*(1-time/10)        */
        else if (model_id == 17)
          {
           dbl Press0, t1, press;
           Press0 = augc[iAC].DataFlt[1];
           t1 = augc[iAC].DataFlt[2];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;

                       if (time_value < 10)
             {
             AC[iAC] = press -Press0*(1.0-exp(-time_value/t1))*(1-time_value/10);
             }
                       else
                         {
             AC[iAC] = press - 0;
                         }
           dAC[iAC][iAC]=1.0;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
             *have_bAC = TRUE;
          }


        /*  Spring condition for microprinting - RBS - 12/12/06 */
        else if (model_id == 20)
          {
           dbl springK[4], forcex, forcey, force, displ;
	   dbl forcey_dcl1, forcey_dcl2, forcex_dcl1, forcex_dcl2;
	   double omega_l, omega_u, t_offset, angle_u;
           int blk_id, ss_id, N_spring;
	   int dcl1_ns, dcl2_ns, kine_ss1, kine_ss2;

	   N_spring = augc[iAC].len_AC - 10;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           for(i=0;i<N_spring;i++){
           	springK[i] = augc[iAC].DataFlt[3+i];
		}
           omega_l  = augc[iAC].DataFlt[3+N_spring];
           omega_u  = augc[iAC].DataFlt[4+N_spring];
           t_offset = augc[iAC].DataFlt[5+N_spring];
           dcl1_ns = ((int) augc[iAC].DataFlt[6+N_spring]);
           dcl2_ns = ((int) augc[iAC].DataFlt[7+N_spring]);
           kine_ss1 = (int) augc[iAC].DataFlt[8+N_spring];
           kine_ss2 = (int) augc[iAC].DataFlt[9+N_spring];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; 
               bAC[iAC][i]=0.0 ; 
		}
               af->Assemble_Jacobian = TRUE;

           forcey = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  4,
                                  "FORCE_Y",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           forcex = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           forcey_dcl1 = evaluate_flux( exo,
                                  dpi,
                                  dcl1_ns,
                                  4,
                                  "FORCE_Y",
                                  blk_id,
                                  kine_ss1,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           forcex_dcl1 = evaluate_flux( exo,
                                  dpi,
                                  dcl1_ns,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  kine_ss1,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           forcey_dcl2 = evaluate_flux( exo,
                                  dpi,
                                  dcl2_ns,
                                  4,
                                  "FORCE_Y",
                                  blk_id,
                                  kine_ss2,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           forcex_dcl2 = evaluate_flux( exo,
                                  dpi,
                                  dcl2_ns,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  kine_ss2,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

	   angle_u = omega_u*(time_value+t_offset) 
			+ 0.5*M_PIE*(1.-omega_u/omega_l);
	   force = (forcey+forcey_dcl1+forcey_dcl2)*sin(angle_u) 
			+ (forcex+forcex_dcl1+forcex_dcl2)*cos(angle_u);
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]= cAC[iAC][i]*sin(angle_u) + bAC[iAC][i]*cos(angle_u); 
               bAC[iAC][i]=0.0 ; 
		}
             *have_cAC = TRUE;

	     *have_dAC = TRUE;

             displ = BC_Types[augc[iAC].BCID].u_BC[augc[iAC].DFID];
           AC[iAC] =  force + springK[0]*displ;
	   if(N_spring >= 2)
		{
           	for(i=1;i<N_spring;i++)
			{
           		AC[iAC] += springK[i]*pow(displ,i+1);
			}
		}
	   dAC[iAC][iAC] = springK[0];
	   if(N_spring >= 2)
		{
           	for(i=1;i<N_spring-1;i++)
			{
	        	dAC[iAC][iAC] += springK[i]*i*pow(displ,i);
			}
		}
          }
        /*  Particle force balance RBS - 5/29/13 */
        else if (model_id == 22)
          {
           dbl forcex, force, force_value, mass;
	   dbl forcex_dcl, angle_u;
           dbl forcex_sens, forcex_dcl_sens;
           int blk_id, ss_id;
	   int dcl_ns, kine_ss;
           int nsp,ns_id,idir,index, index_mesh, k, j;
           dbl xpoint, vpoint, vdotpoint;

           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           dcl_ns = ((int) augc[iAC].DataFlt[3]);
           kine_ss = (int) augc[iAC].DataFlt[4];
           angle_u = augc[iAC].DataFlt[5];
           force_value = augc[iAC].DataFlt[6];
           mass = augc[iAC].DataFlt[7];
           ns_id = 500;
           idir = 1;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; 
               bAC[iAC][i]=0.0 ; 
		}

          af->Assemble_Jacobian = TRUE;
          nsp = match_nsid(ns_id);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
               {
                 k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
                 index = VELOCITY1-1+idir;
                 i = Index_Solution (k, index, 0, 0, blk_id, pg->imtrx);
                 GOMA_EH(i, "Could not resolve index_solution.");
                 vpoint = x[i];
                 vdotpoint = xdot[i];
                 index_mesh = MESH_DISPLACEMENT1-1+idir;
                 i = Index_Solution (k, index_mesh, 0, 0, blk_id, pg->imtrx);
                 GOMA_EH(i, "Could not resolve index_solution.");
                 xpoint = Coor[idir-1][k] + x[i];
                }

           forcex = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

          fprintf(stderr,"forcex done %g %g %g\n",xpoint,vpoint,forcex);

           forcex_dcl = evaluate_flux( exo,
                                  dpi,
                                  dcl_ns,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  kine_ss,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

          fprintf(stderr,"forcex_dcl done %g %g %g\n",xpoint,vpoint,forcex_dcl);

	   force =  (forcex+forcex_dcl);
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]+= bAC[iAC][i]; 
		}
#if 0
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               bAC[iAC][i]=0.0 ; 
		}
           forcex_sens = evaluate_flux_sens( exo,
                                  dpi,
                                  ss_id,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  0,
                                  1,
                                  augc[iAC].BCID,
                                  augc[iAC].DFID,
                                  -1,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);
          fprintf(stderr,"forcex_sens %g %g %g\n",xpoint,vpoint,forcex_sens);
           forcex_dcl_sens = evaluate_flux_sens( exo,
                                  dpi,
                                  dcl_ns,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  kine_ss,
                                  1,
                                  augc[iAC].BCID,
                                  augc[iAC].DFID,
                                  -1,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &bAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

          fprintf(stderr,"forcex_dcl_sens %g %g %g\n",xpoint,vpoint,forcex_dcl_sens);
           dAC[iAC][iAC] = 1*(forcex_sens + forcex_dcl_sens);
	     *have_dAC = TRUE;
/*
             *have_bAC = TRUE;
*/


#endif

             *have_cAC = TRUE;
           AC[iAC] =  force + force_value;
#if 0
fprintf(stderr,"calc %g %g %g %g %g \n",forcex,forcex_dcl,force,forcex_sens,forcex_dcl_sens);
fprintf(stderr,"AC22 %g %g %g %g\n",force,force_value,AC[iAC],BC_Types[9].u_BC[4]);
#endif
          }
/*-------------------------------------------------------------------------*/
        /* Pressure profile                 slah jendoubi                  */
        /* Used to augment pressure boundary conditions                    */
        /* pressure is P=P0*sin(alpha+2*PI*freq*time)                          */
        else if (model_id == 30)
          {
           dbl P0, alpha, freq, press;
           P0 = augc[iAC].DataFlt[1];
           alpha  = augc[iAC].DataFlt[2];
           freq = augc[iAC].DataFlt[3];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;

             AC[iAC] = press - P0*sin(alpha+2*M_PIE*freq*time_value);
           /*P1=P0*sin(alpha+2*M_PIE*freq*time_value);
           printf("time P %lf %lf\n", time_value, P1 );*/

           dAC[iAC][iAC]=1.0;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
             *have_bAC = TRUE;
          }
/*-------------------------------------------------------------------------*/
        /*     Simple constant force (pressure) by  Slah Jendoubi 8/29/2006 */
        /*     Used to augment flow rate boundary conditions                */
        /*     Here the force = F is the force in -x direction              */
        /*     flowrate =0 for F < F0                                       */
        /*     flowrate = a*(F-F0) for F >= F0                              */
        /*     The flowrate is only turned on after time>time1              */
        else if (model_id == 31)
          {
           dbl F, F0, a, force, flowrate, time1;
           int blk_id, ss_id;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           F0 = augc[iAC].DataFlt[3];
           a  = augc[iAC].DataFlt[4];
           time1  = augc[iAC].DataFlt[5];
           flowrate = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

           force = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  3,
                                  "FORCE_X",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);
             F=force;

             if (F < F0)
                 {
                       AC[iAC] = flowrate - 0;
                 }
             else
                 {
                       AC[iAC] = flowrate + a*(F-F0);
                 }

             if (time_value < time1)
                 {
                       AC[iAC] = flowrate - 0;
                 }

           dAC[iAC][iAC]=1.;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
                   *have_bAC = TRUE;

          }
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
        /*     Simple constant force (pressure) by  Slah Jendoubi 8/29/2006 */
        /*     Used to augment flow rate boundary conditions                */
        /*     Here the force = F is the force in -y direction              */
        /*     and p = pressure. to have flow p must be negative            */
        /*     flowrate = 0 for p > p0                                      */
        /*     flowrate = Q0 + a*(p-p0) for p <= p0                         */
        /*     The flowrate is only turned on after time>time1              */
        else if (model_id == 32)
          {
           dbl F, F0, a, p, p0, force, Q, Q0, flowrate, area;
           int blk_id, ss_id;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           p0 = augc[iAC].DataFlt[3];
           a  = augc[iAC].DataFlt[4];
           area  = augc[iAC].DataFlt[5];
           F0 = -area*p0;
           Q0 = 0;
           flowrate = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

           force = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  4,
                                  "FORCE_Y",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           flowrate = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  6,
                                  "VOLUME_FLUX",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);
             F=force;
             p=-F/area;
             Q=flowrate;
          printf("time p VOLUME_FLUX %lf %lf %lf\n", time_value, p, flowrate );

             if (p > p0)
                 {
                       /* do nothing ;*/
                       AC[iAC] = flowrate ;
          /*printf("y32_1 p0 press F VOLUME_FLUX %lf %lf %lf %lf %lf\n", time_value, p0, p, F, Q );*/
                 }
             else
                 {
          /*printf("y32_2 p0 press F VOLUME_FLUX %lf %lf %lf %lf %lf\n", time_value, p0, p, F, Q );*/
                       AC[iAC] = flowrate - (Q0 + a*(p-p0));
                       /*Q = (Q0 + a*(p-p0));*/
          /*printf("time p VOLUME_FLUX %lf %lf %lf\n", time_value, p, Q );*/
                 }
                       dAC[iAC][iAC]=1.;
                       *have_dAC = TRUE;

                       for(i=0;i<NumUnknowns[pg->imtrx];i++){
                               cAC[iAC][i] = 0; }
                         *have_cAC = TRUE;
                       for(i=0;i<NumUnknowns[pg->imtrx];i++){
                               bAC[iAC][i] = 0.0; }
                               *have_bAC = TRUE;
          }
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
        /*     Simple constant Pressure (force) by  Slah Jendoubi 9/06/2006 */
        /*     Used to augment flow rate boundary conditions                */
        /*     Here the force = F is the force in -y direction              */
        /*     Q = flowrate > 0 for all F                                   */
        /*     Q = Q0/F0*F for  Q0 <= Q <= 0                                */
        /*     Q = -a*(F-F0) + Q0 for   Q <= Q0 <= 0                        */
        /*     The flowrate is only turned on after time>time1              */
        /*     F = P1 - P2, P1=0,F>F0>0                                     */
        else if (model_id == 33)
          {
           dbl F0, a, a0, Q, Q0, press, p0, flowrate, time1;
           int blk_id, ss_id;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           F0 = augc[iAC].DataFlt[3];
           a  = augc[iAC].DataFlt[4];
           time1  = augc[iAC].DataFlt[5];
           a0 = a;
           Q0 = a0*F0*0;
           p0 = F0;
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

           flowrate = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  6,
                                  "VOLUME_FLUX",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);
             Q = flowrate;


          /*printf("xtime Q0 F0 press VOLUME_FLUX %lf %lf %lf %lf %lf\n", time_value, Q0, F0, press, Q );*/
             /*if (flowrate > 0)
                 {
                       Q=0;
                 }
             if (flowrate < -1)
                 {
                       Q=-1;
                 }*/

             if (Q > Q0)
                 {
                       /*AC[iAC] = press - Q/a0;*/
                       AC[iAC] = press - ((Q-Q0)/a+p0);
                 }
             else
                 {
                       AC[iAC] = press - ((Q-Q0)/a+p0);
                 }

             if (time_value < time1)
                 {
                 }

           dAC[iAC][iAC]=1.;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
                   *have_bAC = TRUE;

          /*printf("time Q0 F0 press VOLUME_FLUX %lf %lf %lf %lf %lf\n", time_value, Q0, F0, press, Q );*/
          }
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
        /*     Simple constant Pressure (force) by  Slah Jendoubi 9/25/2006 */
        /*     Used to augment flow rate boundary conditions                */
        /*     Here the force = F is the force in -y direction              */
        /*     and the pressure = p.                                        */
        /*                                                                  */
        /*     The folowing is the umbrella constraint (umbrella equation)  */
        /*     Q = flowrate <= 0 for all p                                  */
        /*     Q = flowrate = 0 for all p > p0 (p0 is negavive)             */
        /*     Q = a(p-p0) fot p<=p0                                        */
        /*                                                                  */
        /*     To find the pressure that gives us the right flow rate Q     */
        /*     we use the false-position method                             */
        /*     we assume that p that gives us the right flow rate is bound  */
        /*     between p1 and p2; thus p1 <= p <= p2                        */
        /*     the function that we want to find its zero is f=Q1-Q2        */
        /*     Q1: flow rate given by goma for a fixed pressure             */
        /*     Q2: flow rate given by the umbrella equation                 */
        /*                                                                  */
        else if (model_id == 34)
          {
           static dbl p0, p1, p2, a;
           static dbl p, pl, ph, del, del_p, Q1, Q2, f, fl, fh;
           dbl press, Q, flowrate;
           int blk_id, ss_id;
           static int iteration = 0;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           a  = augc[iAC].DataFlt[3];
           p0 = augc[iAC].DataFlt[4];
           p1 = augc[iAC].DataFlt[5];
           p2 = augc[iAC].DataFlt[6];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
               cAC[iAC][i]=0.0 ; }
               af->Assemble_Jacobian = TRUE;

           flowrate = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  6,
                                  "VOLUME_FLUX",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);
                         Q = flowrate;
             if  (iteration <=6 )
                 {
                         p = p1;
                         pl = p1;

                         Q1 = Q;
                         Q1 = -(p+10);
                         Q2 = a*(p-p0);
                         if (p > p0)
                            {
                               Q2 = 0;
                            }
                         fl = Q1 - Q2;
                         f = fl;
                 }
             else if  (iteration >= 7 && iteration <=9 )
                 {
                         p = p2;
                         ph = p2;
                         Q1 = Q ;
                         Q1 = -(p+10);
                         Q2 = a*(p-p0);
                         if (p > p0)
                            {
                               Q2 = 0;
                            }

                         fh = Q1 - Q2;
                         f = fh;


                         if (fl*fh > 0)
                            {
                               DPRINTF(stderr,"root (pressure) must be bracketed\n");
                               printf("root (pressure) must be bracketed\n");
                            }
                         del_p=ph-pl;
                 }
             else
                 {
                         p = pl + del_p * fl / (fl - fh);
                         Q1 = Q;
                         Q1 = -(p+10);
                         Q2 = a*(p-p0);
                         if (p > p0)
                            {
                               Q2 = 0;
                            }

                         f = Q1 - Q2;

                         if (f < 0)
                            {
                               del = pl - p;
                               pl = p;
                               fl = f;
                            }
                         else
                            {
                               del = ph - p;
                               ph = p;
                               fh = f;
                            }
                         del_p=ph-pl;
                 }
                 printf("iter_p_Q1_Q2_f_Q_fl_fh %d %f %f %f %f %f %f %f\n", iteration, p, Q1, Q2, f, Q, fl, fh);
                 iteration++;
                 AC[iAC] = press - p;


           dAC[iAC][iAC]=1.;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
                   *have_bAC = TRUE;

          }

/*-------------------------------------------------------------------------*/
        /* Pressure profile                 slah jendoubi                  */
        /* Used to augment pressure boundary conditions                    */
        /* pressure is P=P0*sin(alpha+2*PI*freq*time)^2                    */
        else if (model_id == 40)
          {
           dbl P0, alpha, freq, press;
           P0 = augc[iAC].DataFlt[1];
           alpha  = augc[iAC].DataFlt[2];
           freq = augc[iAC].DataFlt[3];
           press = BC_Types[augc[iAC].BCID].BC_Data_Float[augc[iAC].DFID];

              for(i=0;i<NumUnknowns[pg->imtrx];i++){
                       cAC[iAC][i]=0.0 ;         }
                       af->Assemble_Jacobian = TRUE;

             AC[iAC] = press - P0*SQUARE(sin(alpha+2*M_PIE*freq*time_value));
           /*P1=P0*sin(alpha+2*M_PIE*freq*time_value);
           printf("time P %lf %lf\n", time_value, P1 );*/

           dAC[iAC][iAC]=1.0;
           *have_dAC = TRUE;

           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   cAC[iAC][i] = 0; }
             *have_cAC = TRUE;
           for(i=0;i<NumUnknowns[pg->imtrx];i++){
                   bAC[iAC][i] = 0.0; }
             *have_bAC = TRUE;
          }
        /* Simple Nip Force Profile-  RSBay Mar13-06                       */
        /*     Used to augment velocity boundary conditions                */
        /*  Force over area increases linearly from P0 to Pmax in time t1  */
        /*  then decreases back again linearly to time t2 and holds at P0  */
        else if (model_id == 50)
          {
           dbl Press0, NipLoad, PressMax, footprint, web_speed, area, force;
           int blk_id, ss_id;
	   dbl force2, temp[NumUnknowns[pg->imtrx]];
	   int ss_id2;
           blk_id = (int) augc[iAC].DataFlt[1];
           ss_id = (int) augc[iAC].DataFlt[2];
           Press0 = augc[iAC].DataFlt[3];
           NipLoad  = augc[iAC].DataFlt[4];
           footprint = augc[iAC].DataFlt[5];
           web_speed = augc[iAC].DataFlt[6];
           area = augc[iAC].DataFlt[7];

	   memset(cAC[iAC],0, sizeof(double)*NumUnknowns[pg->imtrx]);
	   memset(temp,0, sizeof(double)*NumUnknowns[pg->imtrx]);
           af->Assemble_Jacobian = TRUE;

           force = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  4,
                                  "FORCE_Y",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

           if( augc[iAC].len_AC > 8)
	     {
           ss_id2 = (int) augc[iAC].DataFlt[8];
           force2 = evaluate_flux( exo,
                                  dpi,
                                  ss_id2,
                                  0,
                                  "FORCE_NORMAL",
                                  blk_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  temp,
                                  delta_t,
                                  time_value,
                                  0);
	     }
	   else	{force2 = 0.;}
             *have_cAC = TRUE; 

                       if (time_value < 0)
                         {
             		  target = area* Press0 ;
                         }
                       else if (time_value < footprint/web_speed)
                         {
			  PressMax = 1.5*NipLoad/footprint;
			  target = Press0;
             		  target += PressMax*
				(1. - 4.*SQUARE(web_speed*time_value 
					- 0.5*footprint)/SQUARE(footprint));
			  target *= area;
                         }
                       else
                         {
             		  target = area*Press0;
                         }
           for(i=0;i<NumUnknowns[pg->imtrx];i++){ cAC[iAC][i]  -= temp[i]; } 
             AC[iAC] =  target + force - force2;

          }
        else if(model_id == 51)                 /* Slah user_ac */
        {
          double dist, ratio, ycoord, depth;
          int nsp,ns_id,index, i1, i2;
          ns_id = (int) augc[iAC].DataFlt[1];
          ratio = augc[iAC].DataFlt[2];
          depth = augc[iAC].DataFlt[3];

	  memset(cAC[iAC],0, sizeof(double)*NumUnknowns[pg->imtrx]);  
          af->Assemble_Jacobian = TRUE;

          nsp = match_nsid(ns_id);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
               {
                 k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
                 index = FILL;
                 i1 = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
                 GOMA_EH(i1, "Could not resolve index_solution.");
                 dist = x[i1];
                 index = MESH_DISPLACEMENT2;
                 i2 = Index_Solution (k, index, 0, 0, -1, pg->imtrx);
                 GOMA_EH(i2, "Could not resolve index_solution.");
                 ycoord = Coor[1][k] + x[i2];
                }

          AC[iAC] = ratio*dist - (1.-ratio)*(ycoord+depth);
	  *have_cAC = TRUE;
	  cAC[iAC][i1] = ratio;  
	  cAC[iAC][i2] = -(1.-ratio);  
        }
           /*  radiused upstream bead  */
	else 
	if (model_id == 62)
	{
	  dbl delta_s, radius,dca,sca, lip_angle, alpha1, alpha2;
          dbl xcenter, ycenter, pos_dcl[3]={0,0,0},pos_scl[3]={0,0,0};
          dbl dm_scl[3]={0,0,0};
          int nset_dcl, nset_scl, nsp, i_dcl, i_scl, dir, bc_id;
          nset_dcl = (int) augc[iAC].DataFlt[2];
          nset_scl = (int) augc[iAC].DataFlt[3];
	  dca=M_PIE*augc[iAC].DataFlt[4]/180.0;
	  sca=M_PIE*augc[iAC].DataFlt[5]/180.0;
	  lip_angle=M_PIE*augc[iAC].DataFlt[6]/180.0;
          bc_id = (int) augc[iAC].DataFlt[7];

          alpha1 = 0.5*M_PIE - sca + lip_angle;
          alpha2 = dca - 0.5*M_PIE;
         
          /*  Get DCL and SCL coordinates  */
          nsp = match_nsid(nset_dcl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1, pg->imtrx);
                    GOMA_EH(i, "Could not resolve index_solution.");
                    i_dcl = i;
                    pos_dcl[dir] = Coor[dir][k] + x[i];
                   }
             }
          nsp = match_nsid(nset_scl);
          k = Proc_NS_List[Proc_NS_Pointers[nsp]];
          for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               for (dir = 0; dir < pd->Num_Dim ; dir++)
                   {
                    i = Index_Solution (k, MESH_DISPLACEMENT1+dir, 0, 0, -1, pg->imtrx);
                    GOMA_EH(i, "Could not resolve index_solution.");
                    i_scl = i;
                    pos_scl[dir] = Coor[dir][k];
                    dm_scl[dir] = x[i];
                   }
             }

          /* compute upstream bead radius  */
          radius = ((pos_dcl[0]-pos_scl[0])*sin(-lip_angle)
                          +(pos_dcl[1]-pos_scl[1])*cos(-lip_angle))/
                   (sin(alpha1-lip_angle)-sin(alpha2-lip_angle));
          delta_s = (pos_dcl[0]-pos_scl[0])*cos(-lip_angle)-(pos_dcl[1]-pos_scl[1])
                     *sin(-lip_angle) - radius*(cos(alpha1-lip_angle)-cos(alpha2-lip_angle));

	  xcenter = pos_dcl[0] + radius*cos(alpha2);
	  ycenter = pos_dcl[1] + radius*sin(alpha2);
fprintf(stderr,"AC %g %g %g %g\n",xcenter,ycenter,radius,delta_s);
          BC_Types[bc_id].BC_Data_Float[0] = SQUARE(xcenter);
          BC_Types[bc_id].BC_Data_Float[1] = -2.*(xcenter);
          BC_Types[bc_id].BC_Data_Float[2] = 1.0;
          BC_Types[bc_id+1].BC_Data_Float[0] = SQUARE(ycenter)-SQUARE(radius);
          BC_Types[bc_id+1].BC_Data_Float[1] = -2.*(ycenter);
          BC_Types[bc_id+1].BC_Data_Float[2] = 1.0;

          AC[iAC] = dm_scl[0]+delta_s*cos(lip_angle);
	  *have_dAC = TRUE;
          dAC[iAC][iAC] = cos(lip_angle);
          /*
	  *have_cAC = TRUE;
	  cAC[iAC][i_scl] = 1.0;  
	  cAC[iAC][i_dcl] = cos(-lip_angle);  
	  cAC[iAC][i_dcl] += cos(-lip_angle);  
	  cAC[iAC][i_dcl+1] = -sin(-lip_angle);  
          */
        }  
        else if(model_id == 60) 
        { 
        int nsp,ns_id, num_nodes, ns_id2, mat_id, mn;
        double gap_nom, radius, gap, Vweb, *xpoint, *ypoint;
        double flowrate;
	double p, p3, dns_xpt, cee, dns_ypt, phi[3], phic[3];
        double dxdc, dydc;
        mat_id = (int) augc[iAC].DataFlt[1];
        ns_id = (int) augc[iAC].DataFlt[2];
        ns_id2 = (int) augc[iAC].DataFlt[3];
        mn = map_mat_index(mat_id);
        flowrate = mp_glob[mn]->u_shell_user_par[3];
        gap_nom = mp_glob[mn]->u_shell_user_par[7];
        radius = mp_glob[mn]->u_shell_user_par[5];
        Vweb = mp_glob[mn]->u_shell_user_par[4];
	/*  coordinate of downstream wetting line  */

     	nsp            = match_nsid(ns_id2);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];
        i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
        GOMA_EH(i, "Could not resolve index_solution.");
        dns_xpt = Coor[0][k] + x[i];
        i = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
        GOMA_EH(i, "Could not resolve index_solution.");
        dns_ypt = Coor[1][k] + x[i];

	/*  whole side set list   */
     	nsp            = match_nsid(ns_id);
	num_nodes = Proc_NS_Count[nsp];

	xpoint = smalloc( sizeof(double)*num_nodes);
	ypoint = smalloc( sizeof(double)*num_nodes);
        for (j = 0; j < num_nodes; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint[j] = Coor[0][k] + x[i];
               i = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               ypoint[j] = Coor[1][k] + x[i];
             }
      for ( i = 1; i<num_nodes; i++)
	{
	  p=xpoint[i];
	  p3=ypoint[i];
	  j=i-1;
	  while (j>=0 && (p<xpoint[j]) )
	    {
	      xpoint[j+1]=xpoint[j]; xpoint[j]=p;
	      ypoint[j+1]=ypoint[j]; ypoint[j]=p3;
 	      if( xpoint[i] == xpoint[j]) 
		{
		  fprintf(stderr, "\nMultivalued function detected in list :%d \n",ns_id);
		}
	      j--;
	     }
	 }
    if (num_nodes%2 != 1)
      {
	for( i=num_nodes-1; i >= 0 ; i--)
/*	for( i=0; i < num_nodes ; i++)  */
	{
	  cee = (dns_xpt-xpoint[i])/(xpoint[i+1] - xpoint[i]);
	  if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-1))
	  {
	    phi[0]=-cee+1.;
	    phi[1]=cee;
            phic[0]=-1.;   phic[1]=1.;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1];
            dydc = ypoint[i]*phic[0]+ypoint[i+1]*phic[1];
            dxdc = xpoint[i]*phic[0]+xpoint[i+1]*phic[1];
	    break;
	  }
	}
      }
    else
      {
	for( i=0; i < num_nodes ; i+=2)
        {
          cee = (dns_xpt-xpoint[i])/(xpoint[i+2] - xpoint[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
            phic[0]=4.*cee-3.;  phic[1]=-8.*cee+4.; phic[2]=4.*cee-1.;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1]+ypoint[i+2]*phi[2];
	    dydc = ypoint[i]*phic[0]+ypoint[i+1]*phic[1]+ypoint[i+2]*phic[2];
	    dxdc = xpoint[i]*phic[0]+xpoint[i+1]*phic[1]+xpoint[i+2]*phic[2];
	    break;
	  }
        }
      }
	gap = gap_nom + radius - sqrt(SQUARE(radius)-SQUARE(dns_xpt)) - dns_ypt;
  	AC[iAC] = 0.5*gap*(Vweb*sqrt(1.-SQUARE(dns_xpt)/SQUARE(radius)) +
			Vweb*dxdc/sqrt(SQUARE(dxdc)+SQUARE(dydc))) - flowrate;
        safe_free (xpoint);
        safe_free (ypoint);
	  }
	else if(model_id == 70)
	  {
/* Outlet Condition for roll lamination  - RBS 4/6/2011   */
  	int nsp, ns_id, ns_id2, mat_id, ss_id, unk1, unk2;
        double Vweb;
        double flowrate;
	double dns_ypt2, dns_ypt, dns_xpt;
        mat_id = (int) augc[iAC].DataFlt[1];
        ns_id = (int) augc[iAC].DataFlt[2];
        ns_id2 = (int) augc[iAC].DataFlt[3];
        ss_id = (int) augc[iAC].DataFlt[4];
        Vweb = augc[iAC].DataFlt[5];
	/*  Gap at downstream de-lamination point  */

     	nsp            = match_nsid(ns_id);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];
        unk1 = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
        GOMA_EH(unk1, "Could not resolve index_solution.");
        dns_ypt = Coor[1][k] + x[unk1];
        i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
        GOMA_EH(i, "Could not resolve index_solution.");
        dns_xpt = Coor[0][k] + x[i];
     	nsp            = match_nsid(ns_id2);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];
        unk2 = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
        GOMA_EH(unk2, "Could not resolve index_solution.");
        dns_ypt2= Coor[1][k] + x[unk2];

           flowrate = evaluate_flux( exo,
                                  dpi,
                                  ss_id,
                                  6,
                                  "VOLUME_FLUX",
                                  mat_id,
                                  0,
                                  NULL,
                                  FALSE,
                                  x,
                                  xdot,
                                  &cAC[iAC][0],
                                  delta_t,
                                  time_value,
                                  0);

  	AC[iAC] = flowrate - Vweb*(dns_ypt2 - dns_ypt);
          *have_cAC = TRUE;         
#if 1
 		cAC[iAC][unk1] += Vweb;
 		cAC[iAC][unk2] -= Vweb;
#endif
	  }
	 
	else if (model_id == 80)
	        /* fix a point on the boundary - HvL 10/24/18 */

	{
	int nsp;
        double xpoint;
	
	nsp = match_nsid(103);
      	k = Proc_NS_List[Proc_NS_Pointers[nsp]];
	
	i = Index_Solution(k,MESH_DISPLACEMENT1,0,0,-1, pg->imtrx);
	/*xpoint = Coor[0][k] + x[i];*/

	AC[iAC]=x[i];
	/*AC[iAC]=x[i];*/
	}
	  
/*-------------------------------------------------------------------------*/
    }
  else if (augc[iAC].Type == AC_USERMAT )      /* MT augmenting condition */
    {
        model_id = 0;
        if(augc[iAC].len_AC > 0 )
		{
 		model_id = (int)augc[iAC].DataFlt[0];
		}
	if(model_id == 0)
	  {
  	int nsp,ns_id;
        double pressure,p_value;
        ns_id = (int) augc[iAC].DataFlt[1];
        p_value = augc[iAC].DataFlt[2];
     	nsp            = match_nsid(ns_id);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];

        for (j = 0; j < Proc_NS_Count[nsp]; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, SHELL_LUBP, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               pressure = x[i];
                }
  	   AC[iAC] = pressure - p_value;
	  }
	else if(model_id == 1)
	  {
  	int nsp,ns_id, num_nodes, mn;
        double gap_nom, radius, gap, Vweb, *xpoint, *ypoint;
        double flowrate;
	double p, p3, dns_xpt, cee, dns_ypt, phi[3], phic[3];
        double dxdc, dydc;
        ns_id = (int) augc[iAC].DataFlt[1];
        mn = map_mat_index(augc[iAC].MTID);
        flowrate = mp_glob[mn]->u_shell_user_par[3];
        gap_nom = mp_glob[mn]->u_shell_user_par[7];
        radius = mp_glob[mn]->u_shell_user_par[5];
        Vweb = mp_glob[mn]->u_shell_user_par[4];
     	nsp            = match_nsid(ns_id);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];
	num_nodes = Proc_NS_Count[nsp];

	xpoint = smalloc( sizeof(double)*Proc_NS_Count[nsp]);
	ypoint = smalloc( sizeof(double)*Proc_NS_Count[nsp]);
        for (j = 0; j < num_nodes; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint[j] = Coor[0][k] + x[i];
               i = Index_Solution (k, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               ypoint[j] = Coor[1][k] + x[i];
             }
      for ( i = 1; i<num_nodes; i++)
	{
	  p=xpoint[i];
	  p3=ypoint[i];
	  j=i-1;
	  while (j>=0 && (p<xpoint[j]) )
	    {
	      xpoint[j+1]=xpoint[j]; xpoint[j]=p;
	      ypoint[j+1]=ypoint[j]; ypoint[j]=p3;
 	      if( xpoint[i] == xpoint[j]) 
		{
		  fprintf(stderr, "\nMultivalued function detected in list :%d \n",ns_id);
		}
	      j--;
	     }
	 }
	dns_xpt=augc[iAC].tmp1;
    if (num_nodes%2 != 1)
      {
	for( i=num_nodes-1; i >= 0 ; i--)
/*	for( i=0; i < num_nodes ; i++)  */
	{
	  cee = (dns_xpt-xpoint[i])/(xpoint[i+1] - xpoint[i]);
	  if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-1))
	  {
	    phi[0]=-cee+1.;
	    phi[1]=cee;
            phic[0]=-1.;   phic[1]=1.;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1];
            dydc = ypoint[i]*phic[0]+ypoint[i+1]*phic[1];
            dxdc = xpoint[i]*phic[0]+xpoint[i+1]*phic[1];
	    break;
	  }
	}
      }
    else
      {
	for( i=0; i < num_nodes ; i+=2)
        {
          cee = (dns_xpt-xpoint[i])/(xpoint[i+2] - xpoint[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
            phic[0]=4.*cee-3.;  phic[1]=-8.*cee+4.; phic[2]=4.*cee-1.;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1]+ypoint[i+2]*phi[2];
	    dydc = ypoint[i]*phic[0]+ypoint[i+1]*phic[1]+ypoint[i+2]*phic[2];
	    dxdc = xpoint[i]*phic[0]+xpoint[i+1]*phic[1]+xpoint[i+2]*phic[2];
	    break;
	  }
        }
      }
	gap = gap_nom + radius - sqrt(SQUARE(radius)-SQUARE(dns_xpt)) - dns_ypt;
  	AC[iAC] = 0.5*gap*(Vweb*sqrt(1.-SQUARE(dns_xpt)/SQUARE(radius)) +
			Vweb*dxdc/sqrt(SQUARE(dxdc)+SQUARE(dydc))) - flowrate;
        safe_free (xpoint);
        safe_free (ypoint);
	  }
	else if(model_id == 2)
	  {
  	int nsp,ns_id, num_nodes;
        double p_value, *xpoint, *ypoint;
	double p, p3, dns_xpt, cee, dns_ypt, phi[3];
        ns_id = (int) augc[iAC].DataFlt[1];
        p_value = augc[iAC].DataFlt[2];
     	nsp            = match_nsid(ns_id);
     	k           = Proc_NS_List[Proc_NS_Pointers[nsp]];
	num_nodes = Proc_NS_Count[nsp];

	xpoint = smalloc( sizeof(double)*Proc_NS_Count[nsp]);
	ypoint = smalloc( sizeof(double)*Proc_NS_Count[nsp]);
        for (j = 0; j < num_nodes; j++)
             {
               k = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
               i = Index_Solution (k, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               xpoint[j] = Coor[0][k] + x[i];
               i = Index_Solution (k, SHELL_USER, 0, 0, -1, pg->imtrx);
               GOMA_EH(i, "Could not resolve index_solution.");
               ypoint[j] = x[i];
             }
      for ( i = 1; i<num_nodes; i++)
	{
	  p=xpoint[i];
	  p3=ypoint[i];
	  j=i-1;
	  while (j>=0 && (p<xpoint[j]) )
	    {
	      xpoint[j+1]=xpoint[j]; xpoint[j]=p;
	      ypoint[j+1]=ypoint[j]; ypoint[j]=p3;
 	      if( xpoint[i] == xpoint[j]) 
		{
		  fprintf(stderr, "\nMultivalued function detected in list :%d \n",ns_id);
		}
	      j--;
	     }
	 }
	dns_xpt=augc[iAC].tmp1;
    if (num_nodes%2 != 1)
      {
	for( i=num_nodes-1; i >= 0 ; i--)
/*	for( i=0; i < num_nodes ; i++)  */
	{
	  cee = (dns_xpt-xpoint[i])/(xpoint[i+1] - xpoint[i]);
	  if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-1))
	  {
	    phi[0]=-cee+1.;
	    phi[1]=cee;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1];
	    break;
	  }
	}
      }
    else
      {
	for( i=0; i < num_nodes ; i+=2)
        {
          cee = (dns_xpt-xpoint[i])/(xpoint[i+2] - xpoint[i]);
          if( (cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0)
	      || (cee > 1.0 && i == num_nodes-2))
	  {
	    phi[0]=2.*cee*cee-3.*cee+1.;
	    phi[1]=-4.*cee*cee+4.*cee;
	    phi[2]=2.*cee*cee-cee;
	    dns_ypt = ypoint[i]*phi[0]+ypoint[i+1]*phi[1]+ypoint[i+2]*phi[2];
	    break;
	  }
        }
      }
  	AC[iAC] = dns_ypt - p_value;
	  }
    }

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
                                    const double *const x,
                                    const double *const xdot,
                                    const double delta_t,
                                    const double time_value,
                                    const double *const x_AC,
                                    double *const AC,
                                    double **const cAC,
                                    double **const dAC,
                                    const int numProcUnknowns,
                                    const Exo_DB *const exo,
                                    const Dpi *const dpi,
                                    const Comm_Ex *const cx) {
  /*********************************************************************/
  int i, jAC;
  double inventory, target;
#ifdef PARALLEL
  double global_inventory;
#endif
  inventory = augc[iAC].evol;
#ifdef PARALLEL
  if (Num_Proc > 1) {
    MPI_Allreduce(&inventory, &global_inventory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    inventory = global_inventory;
  }
#endif
  for (i = 0; i < numProcUnknowns; i++) {
    cAC[iAC][i] = augc[iAC].d_evol_dx[i];
  }

  for (jAC = 0; jAC < nAC; jAC++) {
    dAC[iAC][jAC] = 0.0;
  }

  target = augc[iAC].CONSTV;
  AC[iAC] = -target + inventory;
}
/****************************************************************************/
