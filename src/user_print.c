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
 
/* usr_print -- print user-specified information at ea time step to file
 *
 * Typically, this is for pulling nodal point data out into a convenient
 * "x,y" data file format.
 *
 * Created: 1996.02.08 16:39 MST pasacki@sandia.gov
 */
/*
 * "$Id: user_print.c,v 5.1 2007-09-18 18:53:49 prschun Exp $";
 */

#ifdef USE_RCSID
static char rcsid[] =
"$Id: user_print.c,v 5.1 2007-09-18 18:53:49 prschun Exp $";
#endif








#ifdef DEBUG_HKM
#ifdef PARALLEL
#include "mpi.h"
#endif
#endif

#include "usr_print.h"

/*
static int	first_call=TRUE;
static char	sf[] = "%.9e ";	* standard output format for plotting data *
static char	si[] = "%-8d ";	* standard output format for plotting data *
static char	ufn[] = "u_out.d"; * name of user output data file *
FILE	*uf;			* file pointer for this user data *
*/
/*ARGSUSED*/
int
usr_print ( double *t,	            /* time value */
            double dt,              /* time step size */
            double *x,              /* solution vector */
            double **post_proc_vect,
	    int    var)               /* variable of post_proc_vect */
{
  /*  static int status = 0; */
  int retn = 0;

  /*
   * Put a hook in here to print out norms of the solution components
   */
#ifdef DEBUG_HKM
  printf("\tP_%d:  usr_print: Calling Norms routine\n", ProcID);
  fflush(stdout);
#ifdef PARALLEL
  (void) MPI_Barrier(MPI_COMM_WORLD);
#endif
  /* usr_out_hkm(status, *t, dt, x); */
#endif
#ifdef DEBUG_HKM 
#ifdef PARALLEL
  (void) MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  /*
   * Safety catch line -- comment out the line below if you can verify this
   * routine does what you want.
   */
  /* EH(-1, "No usr_print defined."); */

  /*
  if ( first_call )
    {
      if ( (uf = fopen(ufn, "w")) != NULL )
	{
	  first_call = FALSE;
	  fprintf(uf, "# u_out.d");
	}
      else
	{
	  EH(-1, "Could not open user output file.");
	}
    }
    */

  /* 
   * Example:
   *
   * Find the global node number for the desired node that was previously
   * designated as a particular nodeset in the mesh generator.
   *
   *  ns_id          = 2004;
   *  node           = psid2nn(ns_id);
   *  EH(node, "Could not find that nsid.");
   *  initial_pos    = Coor[1][node];
   *  idx            = Index_Solution (node, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
   *  EH(idx, "Could not resolve Index_Solution.");
   *  delta_pos      = x[idx];
   *  actual_pos     = initial_pos + delta_pos;
   *
   *  fprintf(uf, sf, t);
   *  fprintf(uf, sf, dt);
   *  fprintf(uf, sf, actual_pos);
   *  fprintf(uf, "\n");
   *  fflush(uf);			
   */

  /* 
   * Another Example:
   * Print out to file the pairs (arclength, post_proc_vect[node][var]) to
   * a file. This sequence prints out along nodes of an arbitrary node set.
   * This routine should be called at the bottom of rroutine post_process_nodal.
   * WARNING: In this approach the arclength measure can be bad because
   * the nodes on a node set are not necessarily contiguous.  One work
   * around is to output just the x or y coordinate of the node, if one
   * or the other is constant.
   */
/**********Begin example******
  ns_id = 2;
  nsp            = match_nsid(ns_id);  
  initial_pos = 0;
  node           = Proc_NS_List[Proc_NS_Pointers[nsp]];
  fprintf(uf, sf, 0.);
  fprintf(uf, sf, post_proc_vect[var][node]);

  fprintf(uf, "\n");
  fflush(uf);			

  for (j = 1; j < Proc_NS_Count[nsp]; j++)
    {
      node           = Proc_NS_List[Proc_NS_Pointers[nsp]+j];
      node_o         = Proc_NS_List[Proc_NS_Pointers[nsp] + j - 1];
      idx            = Index_Solution (node, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
      idx_o          = Index_Solution (node_o, MESH_DISPLACEMENT1, 0, 0, -1, pg->imtrx);
      EH(idx, "Could not resolve Index_Solution.");
      idy            = Index_Solution (node, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
      idy_o          = Index_Solution (node_o, MESH_DISPLACEMENT2, 0, 0, -1, pg->imtrx);
      EH(idy, "Could not resolve Index_Solution.");
      delta_pos =  sqrt(SQUARE(Coor[0][node]+x[idx]-Coor[0][node_o]-x[idx_o]) 
		      + SQUARE(Coor[1][node]+x[idy]-Coor[1][node_o]-x[idy_o]));
      actual_pos     = initial_pos + delta_pos;
      initial_pos = actual_pos;
      fprintf(uf, sf, actual_pos);
      fprintf(uf, sf, post_proc_vect[var][node]);

      fprintf(uf, "\n");
      fflush(uf);	
    }
*/  /********END EXAMPLE******/		

  /*  return(status); */

  /*  status = 2; */
  return(retn);			/*  Phil's usual default behavior */

} /* END of routine usr_print */
/*****************************************************************************/
/* END of file user_print.c */
/*****************************************************************************/
