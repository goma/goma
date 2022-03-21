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
 *$Id: rf_element_storage.c,v 5.8 2008-10-07 20:32:50 hkmoffa Exp $
 */

//#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "std.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_elem_block_structs.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "rf_allo.h"
#include "rf_element_storage_const.h"
#include "rf_solver.h"
#include "rf_io.h"
#include "mm_eh.h"
#include "mm_input.h"

#include <math.h>

/************************************************************************/
/************************************************************************/
/************************************************************************/

void
setup_element_storage(void)

     /*****************************************************************
      *
      * setup_element_storage()
      *
      *
      *  Loops over element blocks deciding whether to malloc element
      *  storage for that element block.
      *  Note, we could think about different types of element
      *  storage too. However, for now, we just have one type.
      *
      *****************************************************************/
{
  int eb_index, do_malloc, mn;
  ELEM_BLK_STRUCT *eb_ptr; 
  pg->imtrx = 0; // TODO: UCK
  for (eb_index = 0; eb_index < EXO_ptr->num_elem_blocks; eb_index++) {
    do_malloc = FALSE;
    mn = Matilda[eb_index];
    if (mn < 0) {
      continue;
    }
    pd = pd_glob[mn];
    eb_ptr = Element_Blocks + eb_index;
    if (pd->e[pg->imtrx][R_POR_LIQ_PRES]) {
      if (pd->TimeIntegration == TRANSIENT) {
	do_malloc = TRUE;
      }
    }
    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN] || pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      if (pd->TimeIntegration == TRANSIENT) {
	do_malloc = TRUE;
      }
    }
    if(pd->e[pg->imtrx][R_MESH1] && pd->MeshMotion == LAGRANGIAN) {
      if (pd->TimeIntegration == TRANSIENT) {
	do_malloc = TRUE;
      }
    }

    if (do_malloc) {
      init_element_storage(eb_ptr);
      /*Initialize Element-level storage function */
      set_init_Element_Storage(eb_ptr, mn); 

    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
init_element_storage(ELEM_BLK_STRUCT *eb_ptr)

     /*****************************************************************
      *
      * init_element_storage()
      *
      *
      *  Initializes the storage allocation in the Element_Storage
      *  structure. Currently, there is only memory storage for the
      *  the porous flow variable, Sat, here. However, we may expand
      *  the usage of this structure later to include other functionality.
      *  Currently, we malloc this memory in two chunks for space
      *  savings and efficiency. It really doesn't matter as long as 
      *  we know how to unmalloc it.
      *
      *****************************************************************/
{
  int ip_total, i, numStorage;
  ELEMENT_STORAGE_STRUCT *s_ptr;
  double *d_ptr;
  double *base_ptr = NULL;

  if (upd->Total_Num_Matrices > 1) {
    EH(-1, "Element storage not setup to work with multiple matrices.");
  }

  /*
   * Check to make sure that we haven't already allocated storage
   */
  if (!eb_ptr->ElemStorage) {
    /*
     * Determine the number of quadrature points
     */
    ip_total        = eb_ptr->IP_total;
    s_ptr = alloc_struct_1(ELEMENT_STORAGE_STRUCT, 
			   eb_ptr->Num_Elems_In_Block);
    /*
     *  If porous mass lumping is used, we need to store values at
     *  the nodes as well as at the volumetric guass points. We will
     *  do this in a memory inefficient way for now, tacking on the
     *  values at the nodes onto the end of the regular volumetric
     *  quadrature points.
     */
    mp = eb_ptr->MatlProp_ptr;
    if (mp->Porous_Mass_Lump) {
      numStorage = ip_total + eb_ptr->Num_Nodes_Per_Elem;
    } else {
      numStorage = ip_total;
    }

    /*
     * Do a large block allocation for efficiency
     * Argg. See PRS comment in rf_element_storage_struct.h 
     */
     if (pd->e[pg->imtrx][R_POR_LIQ_PRES]) {
       base_ptr = alloc_dbl_1(4 * numStorage * eb_ptr->Num_Elems_In_Block, 
			      DBL_NOINIT); 
     }
     else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN] || pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2] ) {
       base_ptr = alloc_dbl_1(4 * numStorage * eb_ptr->Num_Elems_In_Block, 
			      DBL_NOINIT); 
     }
     else if(pd->e[pg->imtrx][R_MESH1] && pd->MeshMotion == LAGRANGIAN) {

       /* This is for shrinkage stress model for thermexp */
       base_ptr = alloc_dbl_1( numStorage * eb_ptr->Num_Elems_In_Block, 
			      DBL_NOINIT); 
     }
    /*
     * Assign the pointers into the allocated memory block
     */
    d_ptr = base_ptr;
    eb_ptr->ElemStorage = s_ptr;
    for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++, s_ptr++) {
      if (pd->e[pg->imtrx][R_POR_LIQ_PRES] || pd->e[pg->imtrx][R_SHELL_SAT_OPEN] || pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2] ) {
      s_ptr->Sat_QP_tn = d_ptr;
      d_ptr += numStorage;
      s_ptr->p_cap_QP = d_ptr;
      d_ptr += numStorage;
      s_ptr->sat_curve_type = d_ptr;
      d_ptr += numStorage;
      s_ptr->sat_curve_type_old = d_ptr;
      d_ptr += numStorage;
      }
      else if (pd->e[pg->imtrx][R_MESH1] && pd->MeshMotion == LAGRANGIAN) {
        s_ptr->solidified = d_ptr;
	d_ptr += numStorage;
      }
    }
  }

  /* To reference use
   * Element_Blocks->ElemStorage[elem_num].Sat_QP_tn[ip_no]
   */
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
set_init_Element_Storage(ELEM_BLK_STRUCT *eb_ptr, int mn)

     /*****************************************************************
      *
      * set_init_Element_Storage()
      *
      *
      * like its predecessor init_element_storage, this function actually
      * places initial values for the draining and wetting curves for the
      * TANH_HYST function, according to the request in the material property
      * database cards for the current material
      *
      *****************************************************************/
{
  int ip_total, i, j, ifound, ip;
  double sat_switch = 0.0, pc_switch = 0.0, Draining_curve, *ev_tmp;
  int error, num_dim, num_nodes;
  int num_elem, num_elem_blk, num_node_sets, num_side_sets, time_step;
  float	version;		/* version number of EXODUS II */
  int	exoid;			/* ID of the open EXODUS II file */
  char	title[MAX_LINE_LENGTH];	/* title of the EXODUS II database */
  float	ret_float;		/* any returned float */
  char	ret_char[3];		/* any returned character */
  int	num_vars;		/* number of var_type variables */
  char	**var_names = NULL;     /* array containing num_vars variable names */
  char  appended_name[MAX_VAR_NAME_LNGTH];


  /*Quick return if model is not hysteretic in nature */
    
  if(mp_glob[mn]->SaturationModel == TANH_HYST)
    {
      ip_total        = eb_ptr->IP_total;
      Draining_curve   = mp_glob[mn]->u_saturation[8];
      if (Guess_Flag ==4 || Guess_Flag == 5)
	{
	  EH(-1,"Not a smooth restart for hysteretic saturation function. If you really want to do this use read_exoII_file or call us");
	}

      if(Guess_Flag == 5 || Guess_Flag == 6)
	{
	  WH(-1,"Initializing Hysteretic Curve values at all Gauss points with read_exoII_file");
	  CPU_word_size = sizeof(double);
	  IO_word_size  = 0;    

	  exoid = ex_open(ExoAuxFile, EX_READ, &CPU_word_size, &IO_word_size , &version);
	  EH(exoid, "ex_open");

	  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem,
			      &num_elem_blk, &num_node_sets, &num_side_sets);
	  EH(error, "ex_get_init for efv or init guess");

	  /*
	   * Obtain the number of time steps in the exodus file, time_step,
	   * We will read only from the last time step
	   */
	  error = ex_inquire(exoid, EX_INQ_TIME, &time_step, &ret_float, ret_char);
	  EH(error, "ex_inquire");

	  /* Based on problem type and available info in database, extract 
	   * appropriate fields
	   */

	  /*
	   * Get the number of nodal variables in the file, and allocate
	   * space for storage of their names.
	   */
	  error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_vars);
	  EH(error, "ex_get_variable_param");
  
	  /* First extract all nodal variable names in exoII database */
	  if (num_vars > 0) {
	    var_names = alloc_VecFixedStrings(num_vars, (MAX_STR_LENGTH+1));
	    error = ex_get_variable_names(exoid, EX_ELEM_BLOCK, num_vars, var_names);
	    EH(error, "ex_get_variable_names");
	    for (i = 0; i < num_vars; i++) strip(var_names[i]);
	  } else {
	    fprintf(stderr,
		    "Warning: no element variables for saturation stored in exoII input file.\n");
	  }


	  /*****THIS IS WHERE YOU LOAD THEM UP ******/

	  ev_tmp = (double *) smalloc(eb_ptr->Num_Elems_In_Block* sizeof(double));
	  ifound = 0;

	  for(ip = 0; ip < ip_total; ip++)
	    {
	      sprintf(appended_name,  "sat_curve_type%d", ip );
	      for(j=0; j < num_vars; j++)
		{
		  if(!strcasecmp(appended_name,var_names[j]))
		    {
		      /*Found variable so load it into element storage */
		      error  = ex_get_var(exoid, time_step, EX_ELEM_BLOCK, j+1,
					  eb_ptr->Elem_Blk_Id,
					  eb_ptr->Num_Elems_In_Block,
					  ev_tmp);
		      ifound = 1;
		    }
		}
	      if(ifound)
		{
		  for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++) 
		    {
		      eb_ptr->ElemStorage[i].sat_curve_type[ip] = ev_tmp[i];
		    }
		}
	      else
		{
		  EH(-1,"Cannot find an element variable for sat. hysteresis");
		}

	      ifound = 0;

	      sprintf(appended_name,  "sat_switch%d", ip );
	      for(j=0; j < num_vars; j++)
		{
		  if(!strcasecmp(appended_name,var_names[j]))
		    {
		      /*Found variable so load it into element storage */
		      error  = ex_get_var(exoid, time_step, EX_ELEM_BLOCK, j+1,
					  eb_ptr->Elem_Blk_Id,
					  eb_ptr->Num_Elems_In_Block,
					  ev_tmp);
		      ifound = 1;
		    }
		}
	      if(ifound)
		{
		  for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++) 
		    {
		      eb_ptr->ElemStorage[i].Sat_QP_tn[ip] = ev_tmp[i];
		    }
		}
	      else
		{
		  EH(-1,"Cannot find an element variable for sat. hysteresis");
		}

	      ifound = 0;
	      sprintf(appended_name,  "pc_switch%d", ip );
	      for(j=0; j < num_vars; j++)
		{
		  if(!strcasecmp(appended_name,var_names[j]))
		    {
		      /*Found variable so load it into element storage */
		      error  = ex_get_var(exoid, time_step, EX_ELEM_BLOCK, j+1,
					  eb_ptr->Elem_Blk_Id,
					  eb_ptr->Num_Elems_In_Block,
					  ev_tmp);
		      ifound = 1;
		    }
		}
	      if(ifound)
		{
		  for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++) 
		    {
		      eb_ptr->ElemStorage[i].p_cap_QP[ip] = ev_tmp[i];
		    }
		}
	      else
		{
		  EH(-1,"Cannot find an element variable for sat. hysteresis");
		}

	    }
	  
	  error = ex_close(exoid);
	  safer_free((void **) &var_names);
	  free(ev_tmp);
	}
      else /*Initialize as dictated by input cards */
	{

	  if(Draining_curve == 1.0)
	    {
	      sat_switch = mp->u_saturation[0];
	      pc_switch = 1.e-12;
	    }
	  else if (Draining_curve == 0.0)
	    {
	      double sat_max = mp->u_saturation[0];
	      double sat_min = mp->u_saturation[4];
	      double alpha_w = mp->u_saturation[3];
	      double beta_w  = mp->u_saturation[2];

	      pc_switch = 1.e12*alpha_w;
	      sat_switch = sat_max - ( sat_max - sat_min)*0.5*(1.0+tanh( beta_w - alpha_w/pc_switch ) ) ;
	    }
	  else
	    {
	      EH(-1,"TANH_HYST must have 1.0 or 0.0 in  9th spot");
	    }

	  for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++) 
	    {
	      for(ip = 0; ip < ip_total; ip++)
		{
		  eb_ptr->ElemStorage[i].p_cap_QP[ip] = pc_switch;
		  eb_ptr->ElemStorage[i].Sat_QP_tn[ip] = sat_switch;
		  eb_ptr->ElemStorage[i].sat_curve_type[ip] = Draining_curve;
		}
	    }
	}
    }

  if(elc_glob[mn]->thermal_expansion_model == SHRINKAGE)
    {
      ip_total        = eb_ptr->IP_total;
      if (Guess_Flag ==4 || Guess_Flag == 5)
	{
	  EH(-1,"Not a smooth restart for solidification shrinkage model.Use read_exoII_file or call us");
	}

      if(Guess_Flag == 5 || Guess_Flag == 6)
	{
	  EH(-1,"Initializing solidified shrinkage model from exoII file not available yet. Use zero");
	}
	
      // Load em up as all unsolidified

      for(ip = 0; ip < ip_total; ip++)
	{
	  for (i = 0; i < eb_ptr->Num_Elems_In_Block; i++) 
	    {
	      eb_ptr->ElemStorage[i].solidified[ip] = 0.0;
	    }
	}
    }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
free_element_blocks(Exo_DB *exo)

     /*****************************************************************
      *
      * fee_element_blocks()
      *
      *
      *  frees the storage in the Element Blocks structure and all 
      *  underlying storage
      *****************************************************************/
{
  free_element_storage(exo);
  safer_free((void **) &Element_Blocks);      
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
free_element_storage(Exo_DB *exo)

     /*****************************************************************
      *
      * fee_element_storage()
      *
      *
      *  frees the storage in the Element_Storage structures in all
      *  element blocks. It then frees the element blocks themselves
      *****************************************************************/
{
  int eb_index;
  ELEM_BLK_STRUCT *eb_ptr;
  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++) {
    eb_ptr = Element_Blocks + eb_index;
    if (eb_ptr) {
      free_elemStorage(eb_ptr);
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
free_elemStorage(ELEM_BLK_STRUCT *eb_ptr)

     /*****************************************************************
      *
      * free_elemStorage()
      *
      *
      *  frees the storage in the Element_Storage structure.
      *  Note, currently, memory is allocated in two chunks for space
      *  savings. Thus, we free the two chunks here.
      *****************************************************************/
{
  ELEMENT_STORAGE_STRUCT *s_ptr = eb_ptr->ElemStorage;
  if (s_ptr) {
    safer_free((void **) &(s_ptr->Sat_QP_tn));
    safer_free((void **) &(eb_ptr->ElemStorage));
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

double 
get_nodalSat_tnm1_FromES(int lnn)

     /*****************************************************************
      *
      * get_nodalSat_tnm1_FromES()
      *
      *  Gets the t = tnm1 value of the saturation storred in the
      *  element storage at the local node number, lnn.
      *
      *  Note: these accessor functions here and below rely on current
      *        values of  Current_EB_ptr and ei being correct!
      *****************************************************************/
{
  int ielem = ei[pg->imtrx]->ielem;
  int ip_total = Current_EB_ptr->IP_total;
  ELEMENT_STORAGE_STRUCT *es = Current_EB_ptr->ElemStorage + ielem;
  return (es->Sat_QP_tn[ip_total + lnn]);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

double 
get_Sat_tnm1_FromES(int ip)

     /*****************************************************************
      *
      * get_Sat_tnm1_FromES()
      *
      *   Gets the t = tnm1 value of the saturation at the volumetric
      *   quadarture point, ip.
      *
      *****************************************************************/
{
  int ielem = ei[pg->imtrx]->ielem;
  ELEMENT_STORAGE_STRUCT *es = Current_EB_ptr->ElemStorage + ielem;
  return (es->Sat_QP_tn[ip]);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
put_nodalSat_tn_IntoES(int lnn, double sat)

     /*****************************************************************
      *
      * put_nodalSat_tn_IntoES()
      *
      *  Stores the current value of the saturation ( T = tn) into the
      *  element storage at the local node number, lnn.
      *
      *****************************************************************/
{
  int ielem = ei[pg->imtrx]->ielem;
  int ip_total = Current_EB_ptr->IP_total;
  ELEMENT_STORAGE_STRUCT *es = Current_EB_ptr->ElemStorage + ielem;
  es->Sat_QP_tn[ip_total + lnn] = sat;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
put_Sat_tn_IntoES(int ip, double sat)

     /*****************************************************************
      *
      * put_Sat_tn_IntoES()
      *
      *  Stores the current value of the saturation (T = tn) into the
      *  element storage for the current volumetric quadrature point,
      *  ip.
      *
      *****************************************************************/
{
  int ielem = ei[pg->imtrx]->ielem;
  ELEMENT_STORAGE_STRUCT *es = Current_EB_ptr->ElemStorage + ielem;
  es->Sat_QP_tn[ip] = sat;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

