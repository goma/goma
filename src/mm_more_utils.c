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
 *$Id: mm_more_utils.c,v 5.3 2008-03-24 17:43:55 hkmoffa Exp $
 */


/*
 *
 * in this file:
 * 	cnt_nodal_vars()	count TOTAL number of nodal variables for
 *				the problem, including multiplicity of
 *				gas species unknowns. used
 *
 * 	cnt_nodal_vars_post()	count ADDITIONAL number of nodal
 *                              variables to be output as
 *				specified by user in input deck
 *
 *	load_nodal_tkn()	load up information on each of the nodal
 *				variables enumerated above. this information
 *				includes such things as variable names,
 *				variable type, variable kind. Everything is
 *				put into a Results_Description structure.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "std.h"
#include "el_geom.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "el_elm.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "exo_struct.h"		/* defn of Exo_DB */
#include "el_elm_info.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "mm_more_utils.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_bc_const.h"
#include "rf_util.h"

#define GOMA_MM_MORE_UTILS_C

/* cnt_nodal_vars -- count total number of unknown nodal point variables
 *
 * input arguments:
 *	none
 *
 *
 * output arguments:
 * 	none
 *
 * return value:
 *	int		-1	something went wrong in this routine
 *			0,1,... number of nodal point variables
 *
 * Note: We avoid the surface variable unknowns and the pressure for the
 *       moment. Eventually, there should be some mechanism available for
 *	 dumping out the surface species unknown concentrations, as well as
 *       the pressure, that has a different basis function representation
 *	 compared to these conventional unknowns.
 *
 *	 Also, someday it would be nice to have a facility for dumping out
 *	 values of gradients multiplied by thermophysical properties, etc.
 *	 so that vector plots of q=k.grad(T) could be made.
 *
 * Author:	Philip A. Sackinger
 * Created:	Tue Mar 23 07:18:06 MST 1993
 * Revised:	
 */

int 
cnt_nodal_vars(void)
{
  int tnv = 0;
  int v;

#ifdef DEBUG
  static const char yo[] = "cnt_nodal_vars";
#endif

  /* For blot's sake, put the displacements first! */
  for (v = MESH_DISPLACEMENT1; v<(MESH_DISPLACEMENT3+1); v++)
    {
      tnv += goal_post_nodal(v);
    }

  for (v=V_FIRST; v<V_LAST; v++)
    {
      if((v <MESH_DISPLACEMENT1) || (v >MESH_DISPLACEMENT3))
	{
	  tnv += goal_post_nodal(v);
	}
    }

  return(tnv);
}
/* end of cnt_nodal_vars() */

/* cnt_elem_vars -- count total number of unknown element variables
 *
 * input arguments:
 *	none
 *
 *
 * output arguments:
 * 	none
 *
 * return value:
 *	int		-1	something went wrong in this routine
 *			0,1,... number of element variables
 *
 * Author:	Randy Lober
 * Created:	Thurs Aug 13 07:18:06 MST 1998
 * Revised:	
 */

int 
cnt_elem_vars(void)
{
  int   i, j;
  int	tev, *ev_var_mask;
  char	*yo;
  yo     = "cnt_elem_vars";
  tev    = 0;

  ev_var_mask = (int *) smalloc( (V_LAST - V_FIRST)*sizeof(int));
  for ( i = 0; i < V_LAST - V_FIRST; i++ ) {
    ev_var_mask [i] = 0;
  }

  /* Put counter here for array cycling over nodal vars and testing for P0 */
  for (i = 0; i < upd->Num_Mat; i++) {
    for ( j = V_FIRST; j < V_LAST; j++) {
      if ( pd_glob[i]->v[pg->imtrx][j] != V_NOTHING ) {
	if (FALSE && pd_glob[i]->i[pg->imtrx][j] == I_P0) {
	  if (Num_Var_In_Type[pg->imtrx][j] > 1) {
	    fprintf(stderr,
		    "%s: Too many components in variable type for element variable %s (%s)\n",
		    yo,
		    Exo_Var_Names[j].name2,
		    Exo_Var_Names[j].name1 );
	    exit (-1);
	  }
	  if (ev_var_mask[j - V_FIRST] == 0) {
	    /* We just found a candidate for an element variable */
	    tev += Num_Var_In_Type[pg->imtrx][j];
	    ev_var_mask[j - V_FIRST] = 1; /* Only count this variable once */
	  }
        }
	if (FALSE &&pd_glob[i]->i[pg->imtrx][j] == I_P1 ) {	 
	  if (ev_var_mask[j - V_FIRST] == 0) {
	    /* We just found a candidate for an element variable */
	    tev += Num_Var_In_Type[pg->imtrx][j];
	    ev_var_mask[j - V_FIRST] = 1; /* Only count this variable once */
	  }
        }
      }
    }
  }

  safe_free (ev_var_mask);

  return(tev);
}
/* end of cnt_elem_vars() */


/* goal_post_nodal() - code snippet
 *
 *  HKM -> More Documentation would be helpful here
 *
 */

int 
goal_post_nodal(const int var)
{
  int post_flag = FALSE;
  int mat;
  int value;

  for (mat = 0; mat < upd->Num_Mat; mat++)
    {
      post_flag |= ( pd_glob[mat]->i[pg->imtrx][var] == I_Q1 ||
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q2 ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_G ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_G ||
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_GP ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_GP ||
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_GN ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_GN ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_XV ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_XV ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_XG ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_XG ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_HV ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_HG ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_HVG ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_HV ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_HG ||
                     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_HVG ||
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_LSA || 
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q1_D ||
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_D || 
		     pd_glob[mat]->i[pg->imtrx][var] == I_Q2_D_LSA || 
		     pd_glob[mat]->i[pg->imtrx][var] == I_SP);
    }

  if ( post_flag )
  {
    if (var == MASS_FRACTION) {
      value = upd->Max_Num_Species_Eqn;
    } else {
      value = Num_Var_In_Type[pg->imtrx][var];
    }
  }
  else
    {
      value = 0;
    }

  return(value);
}

int 
goal_post_elem(const int var)
{
  int post_flag = FALSE;
  int mat;
  int value;

  char	*yo;
  yo     = "goal_post_elem";

  for (mat = 0; mat < upd->Num_Mat; mat++)
    {
      post_flag |= ( pd_glob[mat]->i[pg->imtrx][var] == I_P0 );
    }

  if ( post_flag )
    {
      if (Num_Var_In_Type[pg->imtrx][var] > 1) {
	fprintf(stderr,
		"%s: Too many components in variable type for element variable\n",
		yo );
	exit (-1);
      }
      value = Num_Var_In_Type[pg->imtrx][var];
    }
  else
    {
      value = 0;
    }

  return(value);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int 
set_nv_tkud(RESULTS_DESCRIPTION_STRUCT *r, /* of just results for Exodus II */
	    const int i,	/* index */
	    const int v,	/* variable index */
	    const int k,	/* kind */
	    const int matIndex, /* Material index for the nodal variable
				 * Pos: Index of the specific material in
				 *      which this variable is defined.
				 * -1 : generic material index
				 * -2 : first unknown at a node having that
				 *      variable type */
	    const char *name,	/* name (short)*/
	    const char *unit,	/* units (unused) */
	    const char *desc,	/* name (long, unused) */
	    const int derivative) /* Does the variable correspond to a time
				     derivative or not */
    /******************************************************************
     *
     * set_nv_tkud()
     *
     * Do frequently-used setup of nodal variables for EXODUS IIv2.02 results
     *
     * return values:
     *
     *		0 == went ok
     *	       -1 == did not
     ******************************************************************/
{
  int status = 0;
  if (i > MAX_NNV) {
    EH(-1, "too many nodal post-process variables for exodus");
  }
  r->nvtype[i] = v;
  r->nvkind[i] = k;
  r->nvmatID[i] = matIndex;
  strcpy(r->nvname[i], name);
  strcpy(r->nvunit[i], unit);
  strcpy(r->nvdesc[i], desc);
  r->nvderivative[i] = derivative;
  return(status);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int 
set_ev_tkud(RESULTS_DESCRIPTION_STRUCT *r, /* elem result of Exodus II       */
	    const int i,	/* index */
	    const int v,	/* variable index */
	    const char *name,	/* name (short)*/
	    const char *unit,	/* units (unused) */
	    const char *desc,	/* name (long, unused) */
	    const int derivative) /* Does the variable correspond to a time
				     derivative or not */    
{
  int status = 0;
  if (i > MAX_NEV) {
    EH(-1, "too many element post-process variables for exodus");
  }
  r->evtype[i] = v;
  strcpy(r->evname[i], name);
  strcpy(r->evunit[i], unit);
  strcpy(r->evdesc[i], desc);
  r->evderivative[i] = derivative;
  return(status);
}

/*
 * Do frequently-used setup of nodal variables for EXODUS IIv2.02 results
 *
 * return values:
 *
 *		0 == went ok
 *	       -1 == did not
 */

int 
load_global_var_info(struct Results_Description *r, /* global results Exodus */
		     const int i, /* index */
		     const char *name) /* name (short)*/
{
  int status;

  status = 0;

  if (r->ngv > MAX_NGV) EH(-1, "too many global post-process varibles for exodus.  Change in rf_solve.c and in names in load_global_var_info");

  strcpy(r->gvname[i], name);

  return(status);
}

void 
sum_total_stress(double sol_vec[],
		 int    var_no,
		 int    k,
		 double nodal_vec[],
		 Exo_DB *exo)

  /*********************************************************************
   *
   * This function gets the nodal stress values for each mode
   * and sums them into a total stress tensor which contains all the
   * mesh nodes, and interpolates the mid-side and centroid values of 
   * all variables with Q1 interpolation on a
   * 9-NODE mesh and interpolates to find their values at the mid-side
   * nodes and at the centroid
   *
   * Now this is set up to be at least compatible w/ parallel computing.
   * We actually load the nodal vector for the current processor only.
   *
   ********************************************************************/
{
  int eb_index;
  int mn, e_start, e_end, ielem, ielem_type, num_local_nodes;
  int iconnect_ptr, var, ktype, i, I, index, ileft, Ileft;

  int iright, Iright;
  int mode;
  
  for (eb_index=0; eb_index<exo->num_elem_blocks; eb_index++)
    {
      mn = Matilda[eb_index];
      pd = pd_glob[mn];
      vn  = vn_glob[mn];

      /* 
       * Assign local pointer pd to appropriate material
       */
      e_start = exo->eb_ptr[eb_index];
      e_end   = exo->eb_ptr[eb_index+1];
      
      for( ielem = e_start; ielem < e_end; ielem++)
	{
	  ielem_type      = Elem_Type(exo, ielem); 

	  num_local_nodes = elem_info(NNODES, ielem_type);
                             /* number of local  basis functions */
    
	  iconnect_ptr    = Proc_Connect_Ptr[ielem];
	                     /* find ptr to beginning of this element's */
     			     /* connectivity list */
	  ktype = k;

	  /* 
	   * First, place the known nodal variable values for this
	   * particular variable and type into the nodal vector.
	   */
	  /* This will zero out the midside node for conjugate problems
	   * at the interface between material.
	   * It needs to be fixed! -RRR
	   */

	  for (i=0; i<num_local_nodes; i++)
	    {
	      I     = Proc_Elem_Connect[iconnect_ptr + i];
	      if(Num_Var_In_Type[pg->imtrx][var_no])
		{
		  index = Index_Solution(I, var_no, ktype, 0, mn, pg->imtrx);
		  if (index != -1)
		    {
		      nodal_vec[I] = sol_vec[index];
		      /* BEWARE: this routine depends upon the exact
		       * ordering in re_fem_const.h. Don't change
		       * it or bad things will happen!
		       */
		      var = var_no + 24;
		      for ( mode=1; mode<vn->modes; mode++)
			{
			  index = Index_Solution(I, var, ktype, 0, mn, pg->imtrx);
			  nodal_vec[I] += sol_vec[index];
			  var += 6;
			}
		    }
		  else
		    {
		      /* Field variable is zero where it is not defined. */
		      nodal_vec[I] = 0.;       
		    }
		}
	    }

	  /*
	   * Rich's famous patch up for lesserly interpolated variables.
	   * Promote quadrilateral Q1 variables to Q2 status, 8 node serendipity
	   * at least, and 9-node biquadratic at best.
	   */

	  /* RRR notes a problem here in 3D.  Should add
	   * check for TRIQUAD_QUAD elem type
	   */

	  if ( pd->v[pg->imtrx][var_no] &&
	       ( pd->i[pg->imtrx][var_no] == I_Q1 ||
                 pd->i[pg->imtrx][var_no] == I_Q1_G ||
		 pd->i[pg->imtrx][var_no] == I_Q1_GP ||
		 pd->i[pg->imtrx][var_no] == I_Q1_GN ||
                 pd->i[pg->imtrx][var_no] == I_Q1_XV ||
                 pd->i[pg->imtrx][var_no] == I_Q1_XG ||
                 pd->i[pg->imtrx][var_no] == I_Q1_HV ||
                 pd->i[pg->imtrx][var_no] == I_Q1_HG ||
                 pd->i[pg->imtrx][var_no] == I_Q1_HVG ||
                 pd->i[pg->imtrx][var_no] == I_SP) && 
	       (ielem_type ==  S_BIQUAD_QUAD || ielem_type == BIQUAD_QUAD ))
	    {
	      /* now interpolate Q1 variables onto 9-node mesh if needed */
	      for (i=4; i<8; i++)
		{
		  I = Proc_Elem_Connect[iconnect_ptr + i]; 
		  /* 
		   * Double check to insure there really are no dofs here.
		   */
		  if (Index_Solution(I,var_no,ktype, 0, mn, pg->imtrx) == -1)
		    {
		      /* 
		       * make node lie halfway between adjacent nodes 
		       * cf. PATRAN local numbering scheme for element.
		       */
		      ileft = i - 4;
		      Ileft = Proc_Elem_Connect[iconnect_ptr + ileft];
		      iright = i - 3;
		      if (iright == 4) iright = 0;
		      Iright = Proc_Elem_Connect[iconnect_ptr + iright];
		      nodal_vec[I] = 
			0.5 * (sol_vec[Index_Solution(Ileft, var_no, ktype, 0, mn, pg->imtrx)] +
			       sol_vec[Index_Solution(Iright, var_no, ktype, 0, mn, pg->imtrx)]);
		      var = var_no + 24;
		      for ( mode=1; mode<vn->modes; mode++)
			{
			  nodal_vec[I] += 
			    0.5 * (sol_vec[Index_Solution(Ileft, var, ktype, 0, mn, pg->imtrx)] +
				   sol_vec[Index_Solution(Iright, var, ktype, 0, mn, pg->imtrx)]);
			  var += 6;
			}
		    }
                  /*
		   * only interpolate centroid in BIQUAD_QUAD
                   */		  
		  if (ielem_type == BIQUAD_QUAD) 
		    {
		      /*
		       *  put centroid in center
		       *  but only if there are no dofs here
		       */
		      I = Proc_Elem_Connect[iconnect_ptr + 8];
        nodal_vec[I] = 0.;
		      if (Index_Solution(I, var_no, ktype, 0, mn, pg->imtrx) == -1)
			{
     for (ileft=0; ileft<4; ileft++)
			    {
			      Ileft = Proc_Elem_Connect[iconnect_ptr + ileft];
			      nodal_vec[I] += 0.25 * 
				sol_vec[Index_Solution(Ileft, var_no, ktype, 0, mn, pg->imtrx)];
			      var = var_no + 24;
			      for ( mode=1; mode<vn->modes; mode++)
				{
				  nodal_vec[I] += 0.25 * 
			 	    sol_vec[Index_Solution(Ileft, var, ktype, 0, mn, pg->imtrx)];
				  var += 6;
				}

			    }
			}
		    }
		}
	    }
	}
    }
  return;
} /* end of sum_total_stress */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
extract_nodal_vec(double sol_vec[], int var_no, int ktype, int matIndex,
		  double nodal_vec[], Exo_DB *exo, int timeDerivative, double time)
    
  /****************************************************************************
   *
   * This function puts the nodal values of the selected variable
   * into a global solution vector which contains all the mesh nodes,
   * and interpolates the mid-side and centroid values of 
   * all variables with Q1 interpolation on a
   * 9-NODE mesh and interpolates to find their values at the mid-side
   * nodes and at the centroid.
   *
   * Now this is set up to be at least compatible w/ parallel computing.
   * We actually load the nodal vector for the current processor only.
   *
   * Written by: Rich Cairncross  27 July 1995
   *
   * Revised: 1997/08/26 14:52 MDT pasacki@sandia.gov
   *
   * Input
   * --------
   *  sol_vec[] = Current global solution vector
   *  var_no    = Variable type to be extracted
   *  ktype     = sub_var number of the variable type to be extracted
   *  matIndex  = material index of the variable to be extracted.
   *             -1 : extract the nonspecific variable
   *             -2 : extract the first variable with var_no no
   *                  matter what material index.
   *  exo       = Exodus database structure
   *  timeDerivative = Are we extracting a time derivative? if so,
   *              then this is true. If not, false.
   *
   * Output
   * -------
   *  nodal_vec[] = nodal vector which receives the value
   *                of the extracted vector. (length number
   *                of nodes)
   *
   *
   *  Special Case Processing:
   *
   *      If var_no is a Species variable, and if we are
   *  solving for nondilute mass transfer, and thus the
   *  number of species may not be equal to the number of species
   *  equations, then this routine can extract the implicit
   *  value of the species variable dropped from the equation
   *  set.
   *    nodal_vec[] in this case for this unknown will be equal
   *  to sum of all of the other species unknowns in the problem.
   *************************************************************************/
{
  int eb_index, mn;
  /*
   * Initialize the output vector to zero. Thus, we don't have to
   * visit each node to get a valid result
   */
  init_vec_value(nodal_vec, 0.0, Num_Node);
  /*
   * Loop backwards through the element blocks
   */
  for (eb_index = exo->num_elem_blocks - 1; eb_index >= 0; eb_index--) {
    mn = Matilda[eb_index];
    pd = pd_glob[mn];
    /*
     * First check to see whether the variable is defined to be
     * in the solution vector in the element block
     * -> If it isn't defined, we don't need to go into the eb routine,
     *    and we will avoid some errors with cross-phase situations.
     */
    if (pd->e[pg->imtrx][var_no]) {
      /*
       *  Next check to see whether we are obtaining material
       *  specific values or general values of the variable
       */
      if (matIndex < 0 || matIndex == mn) {
	/*
	 * Now go into the element block and get the variable at
	 * the nodes.
	 */
	extract_nodal_eb_vec(sol_vec, var_no, ktype, matIndex,
			     eb_index, nodal_vec, exo, timeDerivative, time);
      }
    }
  }
  return;
}
/* end of extract_nodal_vec */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void 
extract_nodal_eb_vec(double sol_vec[], int var_no, int ktype, int matIndex,
		     int eb_index, double nodal_vec[], Exo_DB *exo, 
		     int timeDerivative, double time)

  /******************************************************************************
   *
   * extract_nodal_eb_vec:
   *
   *      This function extracts the a particular nodal variable specified
   *  by the program arguments and loads it up into a global vector.
   *  This routine only operates on one element block.
   *
   * This function puts the nodal values of the selected variable
   * into a global solution vector which contains all the mesh nodes,
   * and interpolates the mid-side and centroid values of 
   * all variables with Q1 interpolation on a
   * 9-NODE mesh and interpolates to find their values at the mid-side
   * nodes and at the centroid
   *
   *  Input
   * --------
   *  sol_vec[] = Current global solution vector
   *  var_no    = Variable type to be extracted
   *  ktype     = sub_var number of the variable type to be extracted
   *  matIndex  = material index of the variable to be extracted.
   *             -1 : extract the nonspecific variable
   *             -2 : extract the first variable with var_no no
   *                  matter what material index.
   *  exo       = Exodus database structure
   *  timeDerivative = Are we extracting a time derivative? if so,
   *              then this is true. If not, false.
   *
   * Output
   * -------
   *  nodal_vec[] = nodal vector which receives the value
   *                of the extracted vector. (length number
   *                of nodes)
   *********************************************************************************/
{
  int mn, e_start, e_end, ielem, ielem_type, num_local_nodes;
  int iconnect_ptr, i, I, index;
  int ileft, Ileft, iright, Iright, midside;
  int lastSpeciesSum, iDof, j;
  MATRL_PROP_STRUCT *matrl;
  double rho;

  /*
   * Find out the material number, mn, from the element block index
   */
  mn = Matilda[eb_index];
  matrl = mp_glob[mn];

  /* 
   * Assign local pointer pd to appropriate material
   */
  pd = pd_glob[mn];

  /*
   *  Check for the existence of a special case for the sum of
   *  the last species constraint condition
   */
  lastSpeciesSum = FALSE;
  if (var_no == MASS_FRACTION) {
    if (matrl->Num_Species_Eqn < matrl->Num_Species) {
      if (ktype == matrl->Num_Species - 1) lastSpeciesSum = TRUE;
    }
  }

  /*
   *  Found out the beginning element number and ending element number
   *  from the element block index
   */
  e_start =  exo->eb_ptr[eb_index];
  e_end   =  exo->eb_ptr[eb_index+1];

  /*
   *  Loop over elements in the element block
   */
  for (ielem = e_start; ielem < e_end; ielem++) {

    /*
     *  Get the element type for this element -> Isn't this the
     *  same for all elements in the element block?
     *  
     */
    ielem_type      = Elem_Type(exo, ielem);

    /*
     * Number of local nodes in the element
     */
    num_local_nodes = elem_info(NNODES, ielem_type);

    /*
     *  find ptr to beginning of this element's connectivity list 
     *
     */
    iconnect_ptr = Proc_Connect_Ptr[ielem]; 

    /* 
     * First, place the known nodal variable values for this
     * particular variable and type into the nodal vector.
     *
     * This will zero out the midside node for conjugate problems
     * at the interface between material.
     * It needs to be fixed! -RRR
     *
     *  The field variable is zero where it is not defined
     */

    for (i = 0; i < num_local_nodes; i++) {
      I = Proc_Elem_Connect[iconnect_ptr + i];
      iDof = 0;
      /*
       * HKM -> Special compatibility section
       */
      if (lastSpeciesSum) {
	index = Index_Solution(I, var_no, 0, iDof, matIndex, pg->imtrx);
	if (index != -1) {
	  nodal_vec[I] = 0.0;
	  for (j = 0; j < matrl->Num_Species_Eqn; j++) {
	    index = Index_Solution(I, var_no, j, iDof, matIndex, pg->imtrx);
	    nodal_vec[I] -= sol_vec[index];
	  }
	  switch (matrl->Species_Var_Type) {
	  case SPECIES_CONCENTRATION:
	      if (matrl->DensityModel == DENSITY_CONSTANT_LAST_CONC) {
		nodal_vec[I] = matrl->u_density[0];
	      } else {
		rho = calc_concentration(matrl, FALSE, NULL);
		nodal_vec[I] += rho;
	      }
	      break;
	  case SPECIES_DENSITY:
	      rho = calc_density(matrl, FALSE, NULL, time);
	      nodal_vec[I] += rho;
	      break;
	  default:
	      if (!timeDerivative) {
		nodal_vec[I] += 1.0;
	      }
	      break;
	  }
	}
      } else {
	index = Index_Solution(I, var_no, ktype, iDof, matIndex, pg->imtrx);
	if (index != -1) {
	  nodal_vec[I] = sol_vec[index];
	} 
      }
    }

    /*
     * Rich's famous patch up for lesserly interpolated variables.
     * Promote quadrilateral Q1 variables to Q2 status, 8 node serendipity
     * at least, and 9-node biquadratic at best.
     */
    /* RRR notes a problem here in 3D. 
     *     Should add check for TRIQUAD_QUAD elem type */

    midside = 0;  
    if (((pd->i[pg->imtrx][var_no] == I_Q1)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_G)   ||
	 (pd->i[pg->imtrx][var_no] == I_Q1_GP)   ||
	 (pd->i[pg->imtrx][var_no] == I_Q1_GN)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_XV)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_XG)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_HV)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_HG)   ||
         (pd->i[pg->imtrx][var_no] == I_Q1_HVG)   ||
	 (pd->i[pg->imtrx][var_no] == I_Q1_D) ||
	 (pd->i[pg->imtrx][var_no] == I_SP)      )
	&& ((ielem_type == S_BIQUAD_QUAD)  ||
	    (ielem_type == BIQUAD_QUAD)         )) {
      midside = 1;
    }
  
    if (midside) {
      /* now interpolate Q1 variables onto 9-node mesh if needed */
      for (i = 4; i < 8; i++) {
	I = Proc_Elem_Connect[iconnect_ptr + i]; 
	/* 
	 * Double check to insure there really are no dofs here.
	 */
	if (Index_Solution(I, var_no, ktype, 0, matIndex, pg->imtrx) == -1) {
	  /* 
	   * make node lie halfway between adjacent nodes 
	   * cf. PATRAN local numbering scheme for element.
	   */
	  ileft = i - 4;
	  Ileft = Proc_Elem_Connect[iconnect_ptr + ileft];
	  iright = i - 3;
	  if (iright == 4) iright = 0;
	  Iright = Proc_Elem_Connect[iconnect_ptr + iright];
	  nodal_vec[I] = 
	      0.5 * (nodal_vec[Ileft] + nodal_vec[Iright]);
#if 0
          if ((pd->i[pg->imtrx][var_no] == I_Q1_HV)   ||
              (pd->i[pg->imtrx][var_no] == I_Q1_HG)   ||
              (pd->i[pg->imtrx][var_no] == I_Q1_HVG)) nodal_vec[I] = -1.;
#endif
	}
      }
      /*
       *  Only interpolate centroid in  BIQUAD_QUAD
       */
      if (ielem_type == BIQUAD_QUAD) {
	I = Proc_Elem_Connect[iconnect_ptr + 8];
	nodal_vec[I] = 0.0; 
	if ( (Index_Solution(I, var_no, ktype, 0, matIndex, pg->imtrx) == -1) ||
             /* for P0 jumps, overwrite jump with interpolant */
             (pd->i[pg->imtrx][var_no] == I_Q1_HV ||
              pd->i[pg->imtrx][var_no] == I_Q1_HG ||
              pd->i[pg->imtrx][var_no] == I_Q1_HVG) ) {
	  for (ileft = 0; ileft < 4; ileft++) {
	    Ileft = Proc_Elem_Connect[iconnect_ptr + ileft];
	    nodal_vec[I] += 0.25 * nodal_vec[Ileft];
	  }
#if 0
          if ((pd->i[pg->imtrx][var_no] == I_Q1_HV)   ||
              (pd->i[pg->imtrx][var_no] == I_Q1_HG)   ||
              (pd->i[pg->imtrx][var_no] == I_Q1_HVG)) nodal_vec[I] = -1.;
#endif
	}     
      }
    } /* END if (midside) */
  } /* END for (ielem = e_start; ielem < e_end; ielem++) */
  return;
}
/********************** end of extract_nodal_eb_vec *****************************/
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void 
extract_elem_vec(const double sol_vec[],
		 const int    ev_indx,
		 const int    var_no,
		 double ***gvec_elem,
		 const Exo_DB *exo )
     
  /***************************************************
   *
   * This function puts the element values of the selected variable
   * into a global solution vector which contains all the elements,
   *
   * Now this is set up to be at least compatible w/ parallel computing.
   * We actually load the nodal vector for the current processor only.
   *
   *
   * Written by: Randy Lober  13 August 1998
   *
   * Revised: 
   *
   ***************************************************/
{
  int eb_index;
  int mn, e_start, e_end, ielem, ielem_type, num_local_nodes;
  int iconnect_ptr, var, ktype, i, I, index;
  int found_quantity;

  for ( eb_index=0; eb_index < exo->num_elem_blocks; eb_index++ ) {

    mn = Matilda[eb_index];
    pd = pd_glob[mn];

    /* 
     * Not all element variables exist in all element blocks.
     * Thus, create_truth_table() didn't malloc all 
     * spots in gvec_elem. For those cases, just skip the 
     * calculation below
     */
    if (gvec_elem[eb_index][ev_indx] == NULL) {
      continue;
    }

    /* 
     * Assign local pointer pd to appropriate material
     */

    e_start = exo->eb_ptr[eb_index];
    e_end   = exo->eb_ptr[eb_index+1];
      
    for (ielem = e_start; ielem < e_end; ielem++) {

      ielem_type      = Elem_Type(exo, ielem); /* func defd in el_geom.h */
      num_local_nodes = elem_info(NNODES, ielem_type);
      iconnect_ptr    = Proc_Connect_Ptr[ielem]; /* find ptr to beginning */
      /* of this element's */
      /* connectivity list */
      var   = var_no;
      ktype = 0;

      /* We're looking at a nodal quantity that should be defined at only 1 
	 node of this element (aka, the pressure value off the hanging center
	 node). If we find it at more than 1 node, we have a serious problem
	 so we're leaving. If we don't find any values, we set the element
	 value to 0. RRl */
      found_quantity = FALSE;
      /* Only do this for elements with a haning center node, otherwise the
	 extraction of the quantity can be ambiguous for non-regular grid models */
      if (ielem_type == BIQUAD_QUAD ||
	  ielem_type == TRIQUAD_HEX ||
	  ielem_type == C_BILINEAR_QUAD ||
	  ielem_type == C_TRILINEAR_HEX ) {
	for (i = 0; i < num_local_nodes; i++) {
	  I     = Proc_Elem_Connect[iconnect_ptr + i];
	  /* NOTE: here, the element variables (such as PRESSURE) are being
	     extracted from the solution vector coming off of the hanging
	     interior nodes, or a given specified node for such a quantity.
	     There should never be more than one of this quantity defined
	     per element, or we have a problem treating it as an element
	     variable. Hence the found_quantity check.                       */
	  index = Index_Solution(I, var, ktype, 0, mn, pg->imtrx);
	  if (index != -1) {
	    /* This should be the one node that has our value - set the element
	       value to this */
	    gvec_elem[eb_index][ev_indx][ielem - e_start] = sol_vec[index];
	    if (found_quantity == TRUE) {
	      fprintf(stderr,
		      "Warning: Too many nodes returning quantities for element variable %s (%s) - may not be accurate\n",
		      Exo_Var_Names[var].name2,
		      Exo_Var_Names[var].name1 );
	      exit (-1);
	    }
	    found_quantity = TRUE;
	  }
	}
      }
      if (found_quantity == FALSE) {
	gvec_elem[eb_index][ev_indx][ielem - e_start] = 0.;   
	/* Field variable is zero where it
	 * is not defined.                 */
      }
#ifdef RRLOBER
      if (found_quantity == FALSE) {
	printf(" No quantity found for variable %s (%s), element %d (%d)\n",
	       Exo_Var_Names[var].name2, Exo_Var_Names[var].name1,
	       ielem, ielem - e_start );
      }
#endif
    }
  }
  return;
}
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/*
 * anneal() - general function maps orignal matl + displacements to new coords
 *
 * Yes, the possibilities are endless with this kind of power...impress me
 * with your imagination...
 *
 * At the least, factor= 1/2 would provide partially relaxed media. Material
 * dependent and coordinate based distortion could provide fancy restart 
 * initial guesses.
 *
 * Created: 1997/08/26 15:18 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void 
anneal_map(const int dim,
	   const double X_old[], 
	   const double displacement[], 
	   double X_new[])
{
  int p;
  double factor=1.0;

  for (p = 0; p < dim; p++) {
      X_new[p] = X_old[p] + factor * displacement[p];
  }
  return;
}
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

int
get_new_coord(double *new_coord[DIM],
	      double *x,
	      const Exo_DB *exo )
{
  int p,i,j;
  int dim = exo->num_dim;
  int num_nodes = exo->num_nodes;
  int displacement_somewhere = FALSE;
  int ln;
  int *moved;
  int e_start = exo->eb_ptr[0];
  int e_end   = exo->eb_ptr[exo->num_elem_blocks];
  int ielem;
  int gnn;
  int var;
  double phi[MDE];

  for(p = 0; p < dim; p++)
    {
      new_coord[p] = (double *) calloc(num_nodes, sizeof(double));
      dcopy1( num_nodes, Coor[p], new_coord[p] );
    }

  for(p = 0; p < upd->Num_Mat; p++)
    displacement_somewhere |= ( pd_glob[p]->e[pg->imtrx][R_MESH1] );

  if ( displacement_somewhere == FALSE ) return (FALSE );

  moved = (int *) calloc( num_nodes, sizeof(int) );


  /*
   * Loop through nodes, find displacement, and add it into 
   * the coordinate
   */


  for(ielem = e_start; ielem < e_end; ielem++)
    {
      double displacement[DIM];

      load_elem_dofptr(ielem, exo, x, x, x, x, 1);

      for(ln = 0; ln < ei[pg->imtrx]->num_local_nodes; ln++)
	{
	  double xi[3] = {0.0, 0.0, 0.0};
	  
	  find_nodal_stu(ln, ei[pg->imtrx]->ielem_type, xi, xi+1, xi+2);

	  gnn = exo->elem_node_list[ exo->elem_node_pntr[ielem] + ln ] ;

	  memset(displacement, 0, sizeof(double)*DIM);
	  
	  if( moved[gnn] != 1 )
	    {
	      for(p = 0; p < DIM; p++)
		{
		  var = MESH_DISPLACEMENT1 + p;
		  
		  for(i = 0; i < ei[pg->imtrx]->dof[var]; i++)
		    {
		      phi[i] = newshape(xi, 
					ei[pg->imtrx]->ielem_type, 
					PSI, 
					ei[pg->imtrx]->dof_list[var][i], 
					ei[pg->imtrx]->ielem_shape,
					pd->i[pg->imtrx][var],
					i);
		    }

		  if( pd->v[pg->imtrx][var] )
		    {
		      for(j = 0; j < ei[pg->imtrx]->dof[var]; j++)
			{
			  displacement[p] += *esp->d[p][j] * phi[j];
			}
		      moved[gnn] = 1;
		    }
		  else
		    displacement[p] = 0.0;
		}
	    }

	  for(p = 0; p < dim; p++)
	    new_coord[p][gnn] += displacement[p];
	}
    }

  safer_free((void **) &moved);
  
  return( TRUE );
}



void
elements_attached_to_NS(int *element_list,
			int NS_ID,
			Exo_DB *exo)
{
	int index;
	int num_nodes;
	int *ns_node_list;
	int i,j,I,e;
	
	if( (index = in_list( NS_ID, 0, exo->ns_node_len, exo->ns_id )) == -1 )
	{	
		EH(-1,"Node set ID not found\n");
	}
	
	num_nodes = exo->ns_num_nodes[index];
	
	ns_node_list = (int *) &( exo->ns_node_list[ exo->ns_node_index[index] ] ) ;				 
					 
	for( i=0, I=0; i< num_nodes; i++)
	{
		I = ns_node_list[i];
		
		for( j=exo->node_elem_pntr[I]; j< exo->node_elem_pntr[I+1]; j++ )
		{
			e = exo->node_elem_list[j];
			
			element_list[e] = TRUE;
		}
	}
	return;
}
/* END of file mm_more_utils.c */
