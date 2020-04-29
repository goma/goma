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
 *$Id: sl_front_setup.c,v 5.3 2009-01-12 16:41:49 hkmoffa Exp $
 */


#include <stdio.h>
#include <string.h>

#include "rf_fem_const.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "sl_rcm.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_node_const.h"
#include "sl_util_structs.h"

extern void setup_front (Exo_DB *);

void
setup_front(Exo_DB *exo)
{
  int i,j,k,l,m, mn, v, ii, kc, idf, nbn_0, nbn, iconnect_ptr, ebi;
  int full_dg, ref_mn=0;
  int index, face, neighbor, m_neighbor, k_neighbor, id_neighbor_centroid, nopdof_offset;
  int one, ccize;
  NODAL_VARS_STRUCT *nv;

  /* Allocate memory for each piece of the fs structure */

  /* Load up required quantities */
  /* viz.
   * fs->bc[ieqn] = 0. for all
   * fs->mdf[I] = Num_Unknowns_Node[I];
   * fs->ncn[ielem] = 1,num_local_nodes ncn += mdf[exo->node_list[ielem]+id]
   * fs->ncod[ieqn] = done in mm_bc.c
   * fs->nop[ielem][id] = Proc_Elem_Connect[ielem] + id
   * fs->nopdof[i*MAX_UNK_ELEM + j] = number of jth dof in element i
   * fs->constraint[i] and fs->el_proc_assign not used right now so just initialized. 
   */

  /*Allocate memory */

  fss->bc    = (double *) array_alloc(1,NumUnknowns[pg->imtrx],sizeof(double));
  fss->ncod  = (int *) array_alloc(1,NumUnknowns[pg->imtrx],sizeof(int));
  fss->ncn   = (int *) array_alloc(1,exo->num_elems,sizeof(int));
  fss->ncn_base = (int *) array_alloc(1,exo->num_elems,sizeof(int));
  fss->mdf = (int *) array_alloc(1,exo->num_nodes,sizeof(int));
  fss->nopp = (int *) array_alloc(1,exo->num_nodes,sizeof(int));
  fss->el_proc_assign = (int *) array_alloc(1,(exo->num_elems + 1),sizeof(int));
  fss->constraint = (int *) array_alloc(1,NumUnknowns[pg->imtrx],sizeof(int));
  fss->level = (int *) array_alloc(1,NumUnknowns[pg->imtrx],sizeof(int));
  fss->perm    = (int *) array_alloc(1, exo->num_elems,sizeof(int));   
  fss->iperm   = (int *) array_alloc(1, exo->num_elems,sizeof(int));   
  fss->mask    = (int *) array_alloc(1, exo->num_elems,sizeof(int));

  /* Some initialization */
  memset(fss->perm, 0, exo->num_elems*sizeof(int)); 
  memset(fss->iperm, 0, exo->num_elems*sizeof(int)); 
  memset(fss->mask, 1,  exo->num_elems*sizeof(int)); 

  /* NB Here the indeces are combined to a single dimensional array with
     a breakout performed in front, due to the difficulty of sending in
     two dimensional arrays into fortran */

  fss->nop = (int *) array_alloc(1,MDE*exo->num_elems, sizeof(int));

  /* For purposes of the global dof connectivity map, you need to test first
     whether to deal with discontinous Galerkin off-element terms for stress*/

  full_dg = 0;
  for (mn = 0; mn < upd->Num_Mat; mn++) 
    {
      if (vn_glob[mn]->dg_J_model == FULL_DG ) 
	{
	  full_dg = 1;
	  ref_mn = mn;
	}
      /*N.B.  Probably ought to record the material in which DGVE is
	active, so you can extract params out of it for the following
	memory allocation and stride.  For now we will do worst case
        scenario and use data off the material ref_mn found here*/
    }

  if(!full_dg)
    {
      fss->nopdof = (int *) array_alloc(1,(MAX_PROB_VAR + MAX_CONC)*MDE*
				             exo->num_elems, sizeof(int));
    }
  else
    {
      /*in the augmented component here the first 4 assumes 4 faces of
	a quadrilateral (2D quads only), the second 4 assumes the number
	of stress components is maximum of 4 for 2D-axi, or 2D.  Need to
	use VIM here.   MDE is the max_number of dofs per element for this
	var.  */
      fss->nopdof = (int *) array_alloc(1,(MAX_PROB_VAR + MAX_CONC)
				       *MDE*exo->num_elems +
                                          4*vn_glob[ref_mn]->modes*4*MDE*
				       exo->num_elems, sizeof(int));
    }

  /* Initialize nop array because some won't get used */

  for (i=0; i<exo->num_elems*MDE; i++)
    {
      fss->nop[i] = 0;
    }
  /*Load up and initialize Obvious quantities first */

  for (i=0; i<NumUnknowns[pg->imtrx]; i++)
    {
      fss->bc[i] = 0.;
      fss->ncod[i] = 0.;
      fss->level[i] = 0;
    }

  for (i = 0; i < exo->num_nodes; i++) {
    nv = Nodes[i]->Nodal_Vars_Info[pg->imtrx];
    fss->mdf[i]=nv->Num_Unknowns;
  }

  for(i=0; i<exo->num_elems; i++) fss->el_proc_assign[i] = 0;
  for(i=0; i<NumUnknowns[pg->imtrx]; i++) fss->constraint[i] = 0;

  fss->nopp[0] = 0;
  for (i=1; i<exo->num_nodes; i++) fss->nopp[i] = fss->nopp[i-1] + fss->mdf[i-1];

  for (i=0; i<exo->num_elems; i++) fss->ncn[i] = fss->ncn_base[i] = 0;


  /* quick test to see if this is all possible */

    for ( ebi=0; ebi<exo->num_elem_blocks-1; ebi++)
    {
      if(elem_info(NNODES,Elem_Type(exo,exo->eb_ptr[ebi])) != 
	 elem_info(NNODES,Elem_Type(exo,exo->eb_ptr[ebi+1])))
	     EH(GOMA_ERROR,"All element blocks must have same number of nodes for this particular frontal solver");
    }


    /* first the traditional  node connectivity map, nop*/

    nbn_0 = elem_info(NNODES,Elem_Type(exo,exo->eb_ptr[0]));   /*big restriction here */
    for( i = 0; i < exo->num_elems; i++) 
      {
	iconnect_ptr = exo->elem_ptr[i];
	for(j=0;j<nbn_0; j++)
	  {
	    fss->nop[nbn_0*i + j] = Proc_Elem_Connect[iconnect_ptr + j];
	  }
      }

    /* Now the global unknown connectivity map */

    /* First do the base, non-discontinuous Galerkin pieces */
    nopdof_offset = 0;
    if(full_dg)nopdof_offset = (4*vn_glob[ref_mn]->modes*4*MDE);

    for( i = 0; i < exo->num_elems; i++) 
      {
	mn = find_mat_number(i, exo);
	  {
	    kc = 0;
	    nbn = elem_info(NNODES,Elem_Type(exo,i));
	    iconnect_ptr = exo->elem_ptr[i];

	    for(j = 0; j < nbn; j++)
	      {
		m = fss->nop[iconnect_ptr + j]; /*current global node num */

		k = fss->nopp[m];
		for(v = V_FIRST; v < V_LAST; v++)
		  {
		    if(v == MASS_FRACTION)
		      {
			idf = Dolphin[pg->imtrx][m][v]*upd->Max_Num_Species_Eqn;
		      }
		    else
		      {
			idf = Dolphin[pg->imtrx][m][v];
		      }
		    if(idf > 0) /*has to be at the node only */
		      /* if(idf > 0 && pd_glob[mn]->v[pg->imtrx][v]) */ /*has to be at the node AND in the block*/
		      {
			for (l=0; l< idf; l++)
			  {
			    ii = k + l;
			    fss->nopdof[i*((MAX_PROB_VAR + MAX_CONC)*MDE +nopdof_offset)+ kc] = ii + 1;
			    kc++;
			    fss->ncn[i]++;
			  }
		      }
		    if(idf > 0) k += idf;
		  }
	      }
	  }
      }

    /* Save base ncn counter for local element, before augmenting */

      for (i=0; i<exo->num_elems; i++) fss->ncn_base[i] = fss->ncn[i];

    /* Now, if necessary add on Discontinous Galerkin pieces */

    for( i = 0; i < exo->num_elems; i++) 
      {
	mn = find_mat_number(i, exo);

	/*bail out if this material doesnt have DG*/  
	/*Actually, this is a problem at interfaces with materials with
	  DG as we are assuming that every element has four connected sides, no? */
	if(vn_glob[mn]->dg_J_model == FULL_DG ) 
	  {
	    kc = 0;
	    nbn = elem_info(NNODES,Elem_Type(exo,i));
	    ei[pg->imtrx]->ielem_type = Elem_Type(exo, i);
	    ei[pg->imtrx]->ielem_shape  = type2shape(ei[pg->imtrx]->ielem_type);
	    ei[pg->imtrx]->num_sides    = shape2sides(ei[pg->imtrx]->ielem_shape);
	    iconnect_ptr = exo->elem_ptr[i];
	    
	    /*unlike above, we deal with the centroid in each element only */
	    j = centroid_node(Elem_Type(exo, i));
	    m = exo->elem_node_list[ exo->elem_node_pntr[i] + j ];
	    
	    k = fss->nopp[m];


	    for (face = 0; face < ei[pg->imtrx]->num_sides; face ++)
	      {
		/*load the neighbor pointer */
		index = exo->elem_elem_pntr[i] + face;
			
		/*load the neighbor element number */
		neighbor = exo->elem_elem_list[index];
		if(neighbor != -1)
		  {  
			
		    /*find the local node number of neighbor corresp. to centroid */
		    id_neighbor_centroid = centroid_node(Elem_Type(exo, neighbor));
		    
		    /*find the global node number of this neighbor */
		    m_neighbor =  exo->elem_node_list[ exo->elem_node_pntr[neighbor] + id_neighbor_centroid ];
		    
		    /*find the global dof number of the first dof at this neighbor centroid */
		    k_neighbor=fss->nopp[m_neighbor];
		    for(v = V_FIRST; v < V_LAST; v++)
		      {
			idf = Dolphin[pg->imtrx][m_neighbor][v];
			if(idf > 0) 
			  {
			    if((v >= POLYMER_STRESS11 && v <= POLYMER_STRESS33) ||
			       (v >= POLYMER_STRESS11_1 && v <= POLYMER_STRESS33_7))
			      { 
				for (l=0; l< idf; l++)
				  {
				    ii = k_neighbor + l;
				    fss->nopdof[i*((MAX_PROB_VAR + MAX_CONC)*MDE 
					          +(4*vn_glob[ref_mn]->modes*4*MDE))+fss->ncn[i]] = ii + 1;
				    fss->ncn[i]++;
				  }
			      }
			    if(idf > 0) k_neighbor += idf;
			  }
		      } 
		  }
	      }  
	  }  /*end if FULL_DG */
	
      }
    
    /* Now some misc. quantities */
    fss->ntra = 0;

    /* Finally, we need an element order map that is bandwidth optimal.  We use
       rcm from aztec, if available. */
      one=1;
      RCM(&one, exo->elem_elem_xadj,
	  exo->elem_elem_adjncy, 
	  fss->mask, fss->perm, &ccize, fss->iperm); 

      printf("rf_solve: used rcm for an optimal order map\n");
      for (i = 0; i < Num_Internal_Elems; i++)  
	{
	  exo->elem_order_map[i]=fss->perm[i];
	}
}


