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
 *  Routines for dumping Goma's Jacobian out to a
 *  machine independent binary file 
 *
 * NOTE:
 *     The  MATRIX_DUMP define must be defined at compile time
 *     for these routines to be compiled.
 */

/*
 *$Id: sl_matrix_dump.c,v 5.2 2009-02-26 23:28:25 hkmoffa Exp $
 */


/* Standard include files */

/*
 *  The xdr functions are in the HPUX_SOURCE namespace.
 */
#ifdef hpux
#ifndef _HPUX_SOURCE
#  define _HPUX_SOURCE
#endif
#endif



#ifdef MATRIX_DUMP
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "mm_as.h"

/*********************** R O U T I N E S  I N   T H I S   F I L E *************
*
*       NAME                            TYPE            CALLED_BY
*
*       matrix_dump_msr ()              void           rf_nonlin_sol
*       matrix_dump_vbr ()              void           rf_nonlin_sol
*
******************************************************************************/



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static unsigned char INT_TO_UCHAR(int value)

/*
*  This function checks whether we can convert to a byte value
*/
  
{
  if (value < 0 || value > 255) {
    fprintf(stderr,
      "matrix_dump ERROR: byte conversion used, but value = %d: Proc = %d\n", 
            value, ProcID);
    exit(-1);
  }
  return((unsigned char) value);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void CHECK_FCLOSE(int result)

/*
*  This function checks the success of a fclose() operations
*/
  
{
  extern int errno;
  if (result) {
    fprintf(stderr,
	    "matrix_dump ERROR: fclose failed Proc = %d: %s\n", ProcID,
	    strerror(errno));
    exit(-1);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void CHECK_XDR(bool result)
  
  /*
  *  This function checks the return status of xdr writes.
  *  We need to be anal about IO checks, because this routine can potentially
  *  put out massive amounts of data, filling up file systems, etc. Also,
  *  there are some MP IO file-locking issues to be handled.
  */
{
  if (!result) {
    fprintf(stderr,
	    "matrix_dump ERROR: xdr write failed Proc = %d\n",
             ProcID);
    exit(-1);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
get_columnNodeList_msr(const double *a, const int *ija, const int ieqn,
		       UMI_LIST_STRUCT *ls_ptr)
{
  int j, ieqnCol, inode = 0;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /*
   * Add the diagonal term
   */
  vd = Index_Solution_Inv(ieqn, &inode, NULL, NULL, NULL, pg->imtrx);
  add_to_umi_int_list(ls_ptr, inode);
  /*
   * Loop through all of the individual terms for the current equation
   * row
   */
  for (j = ija[ieqn]; j < ija[ieqn + 1]; j++) {
    ieqnCol = ija[j];
    vd = Index_Solution_Inv(ieqnCol, &inode, NULL, NULL, NULL, pg->imtrx);
    add_to_umi_int_list(ls_ptr, inode);
  }
  return ls_ptr->Length;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
matrix_dump_msr(struct Aztec_Linear_Solver_System *ams,
		Exo_DB *exo, Dpi *dpi, double *x)

/*****************************************************************************
 *									  
 *  matrix_dump_msr():
 *
 *      This routine will dump a serial file out to disk containing the
 *  Jacobian. The file is meant to be used by the auxilliary program, 
 *  checkJac, to compare to versions of the Jacobian. Ancillary data
 *  meant to enhance the printouts in checkJac are also output to the file.
 * 
 *  The files are named matrix.000, matrix.001, etc.
 *  Overwrites of files are allowed to occur. The files themselves are
 *  written out in native binary format (easy and quick).
 *
 *  The VBR format is used to write files out, even if the internal format
 *  is msr.
 *
 *  The number of matrices to dump out is determined by the value of the
 *  Number_Jac_Dump field in ams. This routine will dump out the first
 *  Number_Jac_Dump matrices constructed by Goma if Number_Jac_Dump is
 *  positive. If Number_Jac_Dump is negative then only the
 *  -Number_Jac_Dump'th jacobian will be dumped out.
 * 
 *****************************************************************************/
{
  extern int  Num_Dim;          /* Number of dimensions in the problem       */
  char   *yo = "matrix_dump_msr ERROR:";
  static int index = 0;         /* Static variable containing the file
                                   index. This routine will dump out files
                                   starting with matrix.000, matrix.001, etc.
                                   into the current directory. Overwriting
                                   of old files is allowed to occur.         */
  char   fname[80];             /* Name of the serial file                   */
  int    file_version = 1;      /* File format version number: if the
				   version number changes, this needs to
				   be bumped. */
  char  *stringPtr;             /* Machine generating this file */
  FILE  *fptr;                  /* File pointer for opening up the serial
                                   file.                                     */
  int *ija;
                                /* Old definitions of VBR index vectors -> 
                                   They are now storred in the AZ_MATRIX     */
  double *a;                    /* Base address of the matrix                */
  int    gNum_Nodes = dpi->num_nodes_global;
                                /* Total global number of nodes in the mesh  */
  int    gNum_Unknowns = dpi->num_dofs_global;
                                /* Global number of unknowns in the Jacobian */
  int    num_owned_nodes = dpi->num_owned_nodes;
                                /* # of unknowns in this mesh updated by 
                                   this processor (internal + border)        */
  int    num_unknowns_node;     /* Number of unknowns at the current
                                   node                                      */
  int    gNodeLocalNext;        /* the global node number of the next
                                   local node to process                     */
  int    gNodeGlobalMin;        /* Global minimum for the global node number
                                   yet to be processed                       */


  int   *gNodeS, *iorder;       /* pointers that will be used to order
                                   block columns in ascending global node
                                   number order                              */

  int   *ordered_gnodes;        /* Pointers that will be used to order nodes
                                   in ascending global node number order     */
  int   *ordered_index;
  NODE_INFO_STRUCT *node;       /* Pointer to the NODE_INFO_STRUCT for the
                                   node corresponding to the current block 
                                   row                                       */
  NODE_INFO_STRUCT *nodeC;      /* Pointer to the NODE_INFO_STRUCT for the
                                   node corresponding to the current block 
                                   column                                    */
  NODE_INFO_STRUCT *nodeT;      /* Temporary pointer                         */
  VARIABLE_DESCRIPTION_STRUCT *vdR, *vdC, *vdT;
  UMI_LIST_STRUCT col_list_umi;
  double **blkMatrix;
  double *a_blk;
  int joffset = 0, ieqnCol;
  char **hndl;

  /*
  *   XDR structures -> allocate space for an XDR structure and assign
  *                     a pointer to it.
  */
  XDR          xdr_structure;
  XDR *xdrs = &xdr_structure;

  /*
  *  A couple of temporary small variables to pack output before writing
  *  to disk
  */
  short int      short_tmp;
  u_short        ushort_tmp;
  float          flt_tmp, flt_coord[3]; 
  unsigned char  byte_tmp;

  int    iblk_row, icol, n1, k, j, i, ieqn_row, num_cols; 
  int    jblk, ieqn, ii, oindex, var_type, sub_type, ib1, itmp;
  int    inode = 0;
  int    i_Var_Desc, i_offset, idof;
  double start_time;
#ifdef DEBUG_HKM
  int    dprint = 0; /* Flag to turn on various amounts of printing */
#endif

  /*
  *   Increment the static index variable. Next time, the matrix
  *   file will have a different name.
  */
  index++;

  /*
   *   Decide on whether to dump out the matrix
   */
  ams->Number_Jac_Dump = 1;
  if (ams->Number_Jac_Dump >= 0) {
    if (ams->Number_Jac_Dump < index) return;
  } else {
    if ((index + ams->Number_Jac_Dump) != 0) return;
  }

  /*
   *      Initialize locations for the VBR index arrays and a[]
   */
  ija   = ams->bindx;
  a     = ams->val;

  /*
  *  Find a name for the matrix dump file  = matrix.000 , matrix.001
  *  Open up the file, and write the header information including
  *  the total global number of nodes onto the file.
  *  If an error occurs, bomb out of the program.
  */
  
  sprintf(fname, "matrix.%03d", index);
#ifdef DEBUG_HKM
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Proc_%d: at top of matrix_dump_msr\n", ProcID); fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  if (ProcID == 0) {
    
    /*
    *  Write an informative message to standard output
    */
    printf("\t\t"); fprint_line(stdout, "=", 64);
    printf("\n\t\tmatrix_dump_msr: dumping jacobian to file %s\n", 
       	   fname);
    start_time = ut();
    
    /*
    *  Write the header information onto the file
    *    - any information not associated with each node in the mesh
    */
    if ((fptr = fopen(fname, "wb")) == NULL) {
      fprintf(stderr, "Couldn't open %s\n", fname);
      exit(-1);
    }
    xdrstdio_create(xdrs, fptr, XDR_ENCODE);
    stringPtr = (char *) array_alloc(1, 80, sizeof(char *));
    strcpy(stringPtr, "goma");
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    CHECK_XDR(xdr_int(xdrs, &file_version));
    (void) gethostname(stringPtr, 80);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    get_date(stringPtr);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    get_time(stringPtr);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    safer_free((void **) &stringPtr);
    CHECK_XDR(xdr_int(xdrs, &gNum_Nodes));
    CHECK_XDR(xdr_int(xdrs, &gNum_Unknowns));
    CHECK_XDR(xdr_int(xdrs, &Num_Dim));
    xdr_destroy(xdrs);
    CHECK_FCLOSE(fclose(fptr));
  }

  /*
  *  Create an index into the Nodes array that sorts the Nodes via their
  *  Global_Node_Num.
  *  In dynamic structures, the nodes are already sorted by their
  *  global node number, so this sort is not needed. ??
  *  In any case, the global node number array, here, has a special entry
  *  added onto the end. The last entry in the array will be global node
  *  number gNumNodes. This will signal the end of processing on this
  *  processor.
  */
  
  ordered_gnodes = (int *) array_alloc(1, 2*num_owned_nodes+2, sizeof(int));
  ordered_index  = ordered_gnodes + num_owned_nodes+1;
  if (ordered_gnodes == NULL) {
    fprintf(stderr, "%s out of memory\n", yo);
    exit(-1);
  }
  
  for (i = 0; i < num_owned_nodes; i++) {
    node = Nodes[i];
    ordered_gnodes[i] = node->Global_Node_Num;
    ordered_index[i]  = i;
  }
  ordered_gnodes[num_owned_nodes] = gNum_Nodes;
  ordered_index[num_owned_nodes]  = num_owned_nodes; 

  /*
   *    Sort the ordered_gnodes[] array in ascending order. The ordered_index
   *  array goes along with the sort. Therefore, it contains an index into
   *  the ordering. (the -1 offset is due to NumRecipes crap)
   */

  sort2_int_int(num_owned_nodes+1, ordered_gnodes-1, ordered_index-1);
#ifdef DEBUG_HKM
  if (ordered_gnodes[num_owned_nodes] != gNum_Nodes) {
    fprintf(stderr,"can't happen\n");
    exit(-1);
  }
  if (ordered_index[num_owned_nodes] != num_owned_nodes) {
    fprintf(stderr,"can't happen\n");
    exit(-1);
  }
#endif

  /*
   * Initialize the umi column list
   */
  col_list_umi.Length = 0;
  col_list_umi.List = NULL;

  /*
  *  Description of the parallel algorithm:
  *
  *  The matrix dump must be sequential in global node number.
  *  Therefore, we will do global gmin_int()'s in order to synchronize 
  *  writing to the same file in MP problems. When a processor has the
  *  lowest global node number remaining, it will process all of its
  *  contiguous global node numbers. The current processor will then
  *  close the file, and participate in the global gmin_int() that
  *  all of the other processors have been waiting on. The signal that
  *  all nodes have been written out is generated when the global 
  *  gmin_int() returns gNum_Nodes.
  */

  /*
  *  Initialize the next global node number counters
  */

  ii = 0;
  oindex = ordered_index[ii];
  gNodeLocalNext = ordered_gnodes[ii];
  gNodeGlobalMin = gmin_int(gNodeLocalNext);

  /*
  *  Loop over the check on whether  gNodeGlobalMin is less than gNum_Nodes
  */

  do {
#ifdef DEBUG_HKM
     if (dprint) print_sync_start(FALSE);
#endif
    /*
    *   Check to see if this processor has the global minimum node number.
    *   If it does, then we need to process the node
    */

    if (gNodeGlobalMin == gNodeLocalNext) {

      /*
      *  Open the output file on the current processor for appending
      */
#ifdef DEBUG_HKM
      if (dprint) {
        printf(
	    "\tProc %d: Writing global nodes starting at %d to file\n",
	    ProcID, gNodeGlobalMin);
      }
#endif
      if ((fptr = fopen(fname, "ab")) == NULL) {
	(void) fprintf(stderr,"Proc %d: Couldn't open %s: %s\n", ProcID, fname,
                       strerror(errno));
	exit(-1);
      }
      xdrstdio_create(xdrs, fptr, XDR_ENCODE);
      
      /*
      *   Loop over nodes -> This processor will try to process as many
      *   nodes as have contiguous global node numbers, now that it owns
      *   the serial file
      */
      
      while ((gNodeGlobalMin == gNodeLocalNext) && 
             (gNodeGlobalMin < gNum_Nodes)) {
	
	/*
	*   Store the pointer to the node structure corresponding to 
	*   the correct global node
	*/
	node = Nodes[oindex];
	
	/*
 	 *   Check to see if this processor actually owns that node
	 */
#ifdef DEBUG_HKM  
	if (! node->Type.Owned) {
	  printf("Error in logic or something\n");
	  exit(-1);
	}
#endif
	/*
	 *    Calculate the block row number of the matrix
	 */
	iblk_row = oindex;
#ifdef DEBUG_HKM
	if (node->Proc_Node_Num != iblk_row) {
	  printf("Screwed up Mesh structure\n");
	  exit(-1);
	}
#endif
	/* 
	*  Value of the beginning row number corresponding to the current
	*  block row. This is also the equation # corresponding to the first
	*  unknown at the node corresponding to the block row.
        *  -> Let's compare against First_Unknown[pg->imtrx] in the node structure
	*     to be sure everything is OK.
	*/
        ieqn_row = node->First_Unknown[pg->imtrx];

	/*
	 *  m1 = number of rows in the current row block This is equal
	 *  to num_unknonws_node, the number of unknowns located at the
	 *  current node.
	 */
	num_unknowns_node = (int) node->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns;	

	/*
	*  Determine the total number of nonzero block columns for the 
	*  current block row and allocate memory for the upcoming sort.
	*/

	for (j = 0, num_cols = 0; j < num_unknowns_node; j++) {
 	  num_cols = get_columnNodeList_msr(a, ija, ieqn_row+j, &col_list_umi);
	}
	
	if (num_cols > 0) {
	  gNodeS = (int *) array_alloc(1, 2*num_cols, sizeof(int));
	  iorder = gNodeS + num_cols;
	  blkMatrix = (double **) array_alloc(1, num_cols, sizeof(double *));
	} else {
	  gNodeS = iorder = NULL;
	  blkMatrix = NULL;
	}
	
	/*
	*  Loop over all the block columns in the current block row. 
	*  Save the global node numbers of the block columns.
	*/
	for (j = 0; j < num_cols; j++) {
	  jblk = col_list_umi.List[j];
	  nodeC = Nodes[jblk];
	  gNodeS[j] = nodeC->Global_Node_Num;
	  iorder[j] = j;
	  n1 = nodeC->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns;
	  blkMatrix[j] = alloc_dbl_1(num_unknowns_node * n1, 0.0);
	}

	/*
	 * Convert the Jacobian entries from msr to vbr for the
	 * current block row.
	 */
	for (i = 0; i < num_unknowns_node; i++) {
          ieqn = ieqn_row + i;

	  /*
	   * Do the diagonal entry first
	   */
	  icol = bin_search_max(col_list_umi.List, num_cols, iblk_row);
	  a_blk = blkMatrix[icol];
	  a_blk[i + i * num_unknowns_node] = a[ieqn];
	  
          /*
	   * Do the off-diagonal entries next
	   */
          inode = iblk_row;
          for (j = ija[ieqn]; j < ija[ieqn + 1]; j++) {
            ieqnCol = ija[j];
            vdC = Index_Solution_Inv(ieqnCol, &inode, NULL, &joffset, NULL, pg->imtrx);
            icol = bin_search_max(col_list_umi.List, num_cols, inode);
	    a_blk = blkMatrix[icol];
	    a_blk[i + joffset * num_unknowns_node] = a[j];
	  }
	}
	
	/*
	 *  Sort the global node numbers of the block columns in increasing
	 *  order
	 */
	sort2_int_int(num_cols, gNodeS-1, iorder-1);
	
	/*
	 *   Print out Header information for the Current Block Row
	 */
#ifdef DEBUG_HKM
	if (dprint > 2) {
	  fprint_line(stdout, "=", 80);
	  printf("Proc: %d Block Row: %d Local Row Number: %d"
		 " Num Rows: %d\n", ProcID, iblk_row, ieqn_row,
		 num_unknowns_node);
	  printf("\t num_block_cols = %d\n", num_cols);
	  fprint_line(stdout,"=", 80);
	}
#endif	
	/*
	*  Dump out Header information for the current block row
	*/
	
	/*       Output the global node number */
        CHECK_XDR(xdr_int(xdrs, &gNodeLocalNext));
	
	/*       Output the processor number */
	CHECK_XDR(xdr_int(xdrs, &ProcID));  
	
	/*       Output the local node number */
	CHECK_XDR(xdr_int(xdrs, &iblk_row));

        /*       Output the coordinates of the current block row
	 *	 node
	 */
        flt_coord[0] = (float) exo->x_coord[iblk_row];
        flt_coord[1] = (float) exo->y_coord[iblk_row];
 	if (Num_Dim > 2) {
          flt_coord[2] = (float) exo->z_coord[iblk_row];
	} else {
	  flt_coord[2] = (float) 0.0;
	}
	CHECK_XDR(xdr_vector(xdrs, (char *) flt_coord,  (u_int) 3, 
                             sizeof(float), (xdrproc_t) xdr_float));  
	
	/*       Output the number of unknowns at this node */
        byte_tmp = INT_TO_UCHAR(num_unknowns_node);
	CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));

	/*
	*    Retrieve the initial equation number for the row and store it 
	*    in ieqn.
	*/
	ieqn = ieqn_row;

        /*
         *   Output the first equation number at the current node
         *  -> assume all equation numbers are contiguous after this one
         */
	CHECK_XDR(xdr_int(xdrs, &ieqn_row));
		
	/*
	 *    Loop over all of the equations at the node, finding the
	 *    variable type, sub_variable type, and checking the node number.
	 */
	for (i = 0; i < num_unknowns_node; i++) {
	  vdT = Index_Solution_Inv(ieqn, &inode, &i_Var_Desc,
		       	           &i_offset, &idof, pg->imtrx);
	  nodeT = Nodes[inode];
#ifdef DEBUG_HKM
	  if (nodeT != node) {
	    fprintf(stderr,"%s Logic error, Proc = %d\n", yo, ProcID);
	    fprintf(stderr,"\t node and nodeT differ %d %d\n",
		    node->Proc_Node_Num, nodeT->Proc_Node_Num);
	    exit (-1);
	  }
#endif
	  /*
	   * Output specific information about the variable
	   */
	  short_tmp = (short int) vdT->Variable_Type;
	  CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	  short_tmp = vdT->MatID;
	  CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	  byte_tmp = INT_TO_UCHAR(vdT->Subvar_Index);
	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
	  byte_tmp = INT_TO_UCHAR(idof);
	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));

	  /*
	   * HKM -> Weights go here. But, will just use a placeholder
	   *        until further research is done
	   */
          flt_tmp = (float) (1.0E-3 * fabs(x[ieqn]) + 1.0E-8);
  	  CHECK_XDR(xdr_float(xdrs, &flt_tmp));
  	  CHECK_XDR(xdr_double(xdrs, x+ieqn));
	  ieqn++;
	}
	
	/*      Output the number of non-zero block columns at this node */
        ushort_tmp = (u_short) num_cols;
        CHECK_XDR(xdr_u_short(xdrs, &ushort_tmp));
	
	/*********************************************************************
         *  Process the block columns in increasing global node number order
	 *********************************************************************/
	
	for (icol = 0; icol < num_cols; icol++) {
	  k = iorder[icol];
	  /*
	   * Find the processor block column (i.e., node number), jblk
	   */
	  jblk = col_list_umi.List[k];

	  /*
	   * Store the location of the block matrix
	   */
          a_blk = blkMatrix[k];
	  
	  /*
	   *  Store the pointer to the node structure for the node corresponding
	   *  to the current block column, nodeC.
	   */
	  
	  nodeC = Nodes[jblk]; 
#ifdef DEBUG_HKM
	  if (nodeC->Global_Node_Num != gNodeS[icol]) {
	    fprintf(stderr,"%s Logic error\n", yo);
	    fprintf(stderr,
		    "nodeC->Global_Node_Num and gNodeS[icol] differ: %d %d \n",
		    nodeC->Global_Node_Num, gNodeS[icol]);
	    exit(-1);
	  }
#endif  
	  /*
	   *  Calculate the equation number for the first unknown 
	   *  in the current block column
	   */
          ib1 = nodeC->First_Unknown[pg->imtrx];
	  
	  /*
	   *   Calculate the number of columns in the current block column
	   *   and cross check.
	   */
	  n1 = nodeC->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns;

	  /*
	   *     Print out an ascii header information
	   *     -> Mostly for debugging purposes.
	   */
#ifdef DEBUG_HKM 
	  if (dprint > 3) {
	    printf("\t Block Row: %d Block Column: %d"
		   " Num Cols: %d\n", iblk_row, jblk, n1);
	  }
	  if (dprint > 5) {
	    fprint_line(stdout, "-", 80);
	    for (i = 0; i < num_unknowns_node; i++) {
	      for (j = 0; j < n1; j++)
		  printf("a[%d]: %e ", j*num_unknowns_node+i,
			 a_blk[j*num_unknowns_node+i]);
	      printf("\n");
	    }
	    fprint_line(stdout, "-", 80);
	  }
#endif
	  /*
	   *  Dump out the header information for the node corresponding to the
           *  block column.
	   */
	  
	  /*        Output global node number of the block column */
  	  CHECK_XDR(xdr_int(xdrs, &gNodeS[icol]));
	  
	  /*        Output the local node number of the block column */
  	  CHECK_XDR(xdr_int(xdrs, &jblk));

          /* 
           *      Output the coordinates of the node corresponding to the
           *      current column 
           */
          flt_coord[0] = (float) exo->x_coord[jblk];
          flt_coord[1] = (float) exo->y_coord[jblk];
	  if (Num_Dim > 2) {
            flt_coord[2] = (float) exo->z_coord[jblk];
	  } else {
	    flt_coord[2] = 0.0;
	  }
          CHECK_XDR(xdr_vector(xdrs, (char *) flt_coord, (u_int) 3, 
                               sizeof(float), (xdrproc_t) xdr_float));  
	  
	  /*        Output the number of unknowns at the column global node */
          byte_tmp = INT_TO_UCHAR(n1);
  	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
	  
	  /*
	  *         Output the identity of first unknown at the current
          *         block column
	  */
          ieqn = nodeC->First_Unknown[pg->imtrx];
  	  CHECK_XDR(xdr_int(xdrs, &ieqn));

	  /*
	  *         Output the identity of the unknowns for the current
	  *         block column
	  */
	  for (i = 0; i < n1; i++) {
            vdT = Index_Solution_Inv(ieqn, &inode, &i_Var_Desc,
		                     &i_offset, &idof, pg->imtrx);
#ifdef DEBUG_HKM
	    nodeT = Nodes[inode];
	    if (nodeT != nodeC) {
	      fprintf(stderr,"%s Logic error, Proc = %d\n", yo, ProcID);
	      fprintf(stderr,"\t nodeC and nodeT differ %d %d\n",
		      nodeC->Proc_Node_Num, nodeT->Proc_Node_Num);
	      exit (-1);
	    }
#endif
	    /*
	     * Output specific information about the variable
	     */
	    short_tmp = (short int) vdT->Variable_Type;
	    CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	    short_tmp = vdT->MatID;
  	    CHECK_XDR(xdr_short(xdrs,  &short_tmp));
            byte_tmp = INT_TO_UCHAR(vdT->Subvar_Index);
  	    CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
            byte_tmp = INT_TO_UCHAR(idof);
  	    CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
    
	    ieqn++;
	  }
	  
	  /*
	   *  Dump out the matrix for the current block-row block-column
	   */
	  itmp = num_unknowns_node*n1;
          CHECK_XDR(xdr_vector(xdrs, (char *) a_blk, (u_int) itmp,
                               sizeof(double), (xdrproc_t) xdr_double));
	  
	} /* END for (icol = 0; icol < num_cols - loop over block cols */
	/*********************************************************************/
	/*
	*   Free allocated memory, allocated within the loop over global
	*   node numbers.
	*/

	for (j = 0; j < num_cols; j++) {
          safer_free((void **) (blkMatrix + j));
	}
	safer_free((void **) &blkMatrix);
        free_umi_list(&col_list_umi);
	safer_free((void **) &gNodeS);

	
        /*
	 *  Increment the node on this processor to the node with the
	 *  next lowest global node number. End conditions are already
         *  taken care of.
	 */
	ii++;
	oindex = ordered_index[ii];
	gNodeLocalNext = ordered_gnodes[ii];

        /*
         *  Increment the global minimum node number
         */
        gNodeGlobalMin++;
	
      } /* End of while loop over processing the local node numbers */

      /*
       *  Close the output file on that processor, and thereby flush
       *  any buffered output
       */
#ifdef DEBUG_HKM
     if (dprint) {
       printf(
 "\tProc = %d, Giving up control: gNodeLocalNext = %d, gNodeGlobalMin = %d\n",
              ProcID, gNodeLocalNext, gNodeGlobalMin);
     }
#endif
      xdr_destroy(xdrs);
      CHECK_FCLOSE(fclose(fptr));
      
    } /* END of -> if (gNodeGlobalMin == gNodeLocalNext) { <- */
    
    /*
    *  Go get the next global minimum node number. And, check to see
    *  whether there still is a processor with more nodes to write out.
    *  If there is, go to the top of the routine and start again.
    */

#ifdef DEBUG_HKM
    if (dprint) {
      if (dprint > 1)
	printf("\t\tProc = %d, gNodeLocalNext = %d\n", ProcID, 
	       gNodeLocalNext);
      print_sync_end(FALSE);
    }
#endif     
    gNodeGlobalMin = gmin_int(gNodeLocalNext);
  } while(gNodeGlobalMin <gNum_Nodes);
  
  /*
   *  Free memory allocated at the beginning of the routine
   */
  safer_free((void **) &ordered_gnodes);
  
  /*
   *  Write an informative message to standard output
   */
  if (ProcID == 0) {
    start_time = ut() - start_time;
    printf("\n\t\tmatrix_dump_msr: finished, dump time =  %g sec\n", 
	   start_time);
    printf("\t\t"); fprint_line(stdout, "=",64);
  }
} /********************** End matrix_dump_msr *******************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
matrix_dump_vbr(struct Aztec_Linear_Solver_System *ams,
		Exo_DB *exo, Dpi *dpi, double *x)

/*****************************************************************************
 *									  
 *  matrix_dump_vbr():
 *
 *      This routine will dump a serial file out to disk containing the
 *  Jacobian. The file is meant to be used by the auxilliary program, 
 *  checkJac, to compare to versions of the Jacobian. Ancillary data
 *  meant to enhance the printouts in checkJac are also output to the file.
 * 
 *  The files are named matrix.000, matrix.001, etc.
 *  Overwrites of files are allowed to occur. The files themselves are
 *  written out in native binary format (easy and quick).
 *
 *  An annotated VBR format is used to write files out.
 *
 *  The number of matrices to dump out is determined by the value of the
 *  Number_Jac_Dump field in ams. This routine will dump out the first
 *  Number_Jac_Dump matrices constructed by Goma if Number_Jac_Dump is
 *  positive. If Number_Jac_Dump is negative then only the
 *  -Number_Jac_Dump'th jacobian will be dumped out.
 * 
 *****************************************************************************/
{
  extern int  Num_Dim;          /* Number of dimensions in the problem       */
  char   yo[] = "matrix_dump_vbr ERROR:";
  static int index = 0;         /* Static variable containing the file
                                   index. This routine will dump out files
                                   starting with matrix.000, matrix.001, etc.
                                   into the current directory. Overwriting
                                   of old files is allowed to occur.         */
  char   fname[80];             /* Name of the serial file                   */
  char  *stringPtr;
  int    file_version = 1;
  FILE  *fptr;                  /* File pointer for opening up the serial
                                   file.                                     */
  int *indx, *bindx, *rpntr, *bpntr, *cpntr;
                                /* Old definitions of VBR index vectors -> 
                                   They are now storred in the AZ_MATRIX     */
  double *a;                    /* Base address of the matrix                */
  double *a_tmp;                /* Temporary pointer into a                  */
  int    gNum_Nodes = dpi->num_nodes_global;
                                /* Total global number of nodes in the mesh  */
  int    gNum_Unknowns = dpi->num_dofs_global;
                                /* Global number of unknowns in the Jacobian */
  int    num_owned_nodes = dpi->num_owned_nodes;
                                /* # of unknowns in this mesh updated by 
                                   this processor (internal + border)        */
  int    num_unknowns_node;     /* Number of unknowns at the current
                                   node                                      */
  int    gNodeLocalNext;        /* the global node number of the next
                                   local node to process                     */
  int    gNodeGlobalMin;        /* Global minimum for the global node number
                                   yet to be processed                       */


  int   *gNodeS, *iorder;       /* pointers that will be used to order
                                   block columns in ascending global node
                                   number order                              */

  int   *ordered_gnodes;        /* Pointers that will be used to order nodes
                                   in ascending global node number order     */
  int   *ordered_index;
  NODE_INFO_STRUCT *node;       /* Pointer to the NODE_INFO_STRUCT for the
                                   node corresponding to the current block 
                                   row                                       */
  NODE_INFO_STRUCT *nodeC;      /* Pointer to the NODE_INFO_STRUCT for the
                                   node corresponding to the current block 
                                   column                                    */
  NODE_INFO_STRUCT *nodeT;      /* Temporary pointer                         */
  VARIABLE_DESCRIPTION_STRUCT *vdR, *vdC, *vdT;

  /*
  *   XDR structures -> allocate space for an XDR structure and assign
  *                     a pointer to it.
  */
  XDR          xdr_structure;
  XDR *xdrs = &xdr_structure;

  /*
  *  A couple of temporary small variables to pack output before writing
  *  to disk
  */
  short int      short_tmp;
  u_short        ushort_tmp;
  float          flt_tmp, flt_coord[3]; 
  unsigned char  byte_tmp;

  int    iblk_row, icol, n1, m1,  k, j, i, ieqn_row, ival, num_cols; 
  int    jblk, ieqn, ii, oindex, var_type, sub_type, ib1, ib2, itmp;
  int    inode = -1, i_Var_Desc = -1, i_offset = -1, idof = -1;
  double start_time;
#ifdef DEBUG_HKM
  int    dprint = 0; /* Flag to turn on various amounts of printing */
#endif

  /*
  *   Increment the static index variable. Next time, the matrix
  *   file will have a different name.
  */
  index++;

  /*
   *   Decide on whether to dump out the matrix
   */
  if (ams->Number_Jac_Dump >= 0) {
    if (ams->Number_Jac_Dump < index) return;
  } else {
    if ((index + ams->Number_Jac_Dump) != 0) return;
  }

  /*
   *      Initialize locations for the VBR index arrays and a[]
   */
  rpntr = ams->rpntr;
  indx  = ams->indx;
  bpntr = ams->bpntr;
  bindx = ams->bindx;
  cpntr = ams->cpntr;
  a     = ams->val;

  /*
  *  Find a name for the matrix dump file  = matrix.000 , matrix.001
  *  Open up the file, and write the header information including
  *  the total global number of nodes onto the file.
  *  If an error occurs, bomb out of the program.
  */
  
  sprintf(fname, "matrix.%03d", index);
  if (ProcID == 0) {
    
    /*
    *  Write an informative message to standard output
    */
    printf("\t\t"); fprint_line(stdout, "=",64);
    printf("\n\t\tmatrix_dump_vbr: dumping jacobian to file %s\n", 
       	   fname);
    start_time = ut();
    
    /*
    *  Write the header information onto the file
    *    - any information not associated with each node in the mesh
    */
    if ((fptr = fopen (fname, "wb")) == NULL) {
      fprintf (stderr, "Couldn't open %s\n", fname);
      exit (-1);
    }

    xdrstdio_create(xdrs, fptr, XDR_ENCODE);
    stringPtr = (char *) array_alloc(1, 80, sizeof(char *));
    strcpy(stringPtr, "goma (vbr)");
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    CHECK_XDR(xdr_int(xdrs, &file_version));
    (void) gethostname(stringPtr, 80);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    get_date(stringPtr);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    get_time(stringPtr);
    CHECK_XDR(xdr_wrapstring(xdrs, &stringPtr));
    safer_free((void **) &stringPtr);
    CHECK_XDR(xdr_int(xdrs, &gNum_Nodes));
    CHECK_XDR(xdr_int(xdrs, &gNum_Unknowns));
    CHECK_XDR(xdr_int(xdrs, &Num_Dim));
    xdr_destroy(xdrs);
    CHECK_FCLOSE(fclose(fptr));
  }

  /*
  *  Create an index into the Nodes array that sorts the Nodes via their
  *  Global_Node_Num.
  *  In dynamic structures, the nodes are already sorted by their
  *  global node number, so this sort is not needed. ??
  *  In any case, the global node number array, here, has a special entry
  *  added onto the end. The last entry in the array will be global node
  *  number gNumNodes. This will signal the end of processing on this
  *  processor.
  */
  
  ordered_gnodes = (int *) array_alloc(1, 2*num_owned_nodes+2, sizeof(int));
  ordered_index  = ordered_gnodes + num_owned_nodes+1;
  if (ordered_gnodes == NULL) {
    fprintf(stderr, "%s out of memory\n", yo);
    exit(-1);
  }
  
  for (i = 0; i < num_owned_nodes; i++) {
    node = Nodes[i];
    ordered_gnodes[i] = node->Global_Node_Num;
    ordered_index[i]  = i;
  }
  ordered_gnodes[num_owned_nodes] = gNum_Nodes;
  ordered_index[num_owned_nodes]  = num_owned_nodes; 

  /*
  *    Sort the ordered_gnodes[] array in ascending order. The ordered_index
  *  array goes along with the sort. Therefore, it contains an index into
  *  the ordering. (the -1 offset is due to NumRecipes crap)
  */

  sort2_int_int(num_owned_nodes+1, ordered_gnodes-1, ordered_index-1);
#ifdef DEBUG_HKM
  if (ordered_gnodes[num_owned_nodes] != gNum_Nodes) {
    fprintf(stderr,"can't happen\n");
    exit(-1);
  }
  if (ordered_index[num_owned_nodes] != num_owned_nodes) {
    fprintf(stderr,"can't happen\n");
    exit(-1);
  }
#endif

  /*
  *  Description of the parallel algorithm:
  *
  *  The matrix dump must be sequential in global node number.
  *  Therefore, we will do global gmin_int()'s in order to synchronize 
  *  writing to the same file in MP problems. When a processor has the
  *  lowest global node number remaining, it will process all of its
  *  contiguous global node numbers. The current processor will then
  *  close the file, and participate in the global gmin_int() that
  *  all of the other processors have been waiting on. The signal that
  *  all nodes have been written out is generated when the global 
  *  gmin_int() returns gNum_Nodes.
  */

  /*
  *  Initialize the next global node number counters
  */

  ii = 0;
  oindex = ordered_index[ii];
  gNodeLocalNext = ordered_gnodes[ii];
  gNodeGlobalMin = gmin_int(gNodeLocalNext);

  /*
  *  Loop over the check on whether  gNodeGlobalMin is less than gNum_Nodes
  */

  do {
#ifdef DEBUG_HKM
     if (dprint) print_sync_start(FALSE);
#endif
    /*
    *   Check to see if this processor has the global minimum node number.
    *   If it does, then we need to process the node
    */

    if (gNodeGlobalMin == gNodeLocalNext) {

      /*
      *  Open the output file on the current processor for appending
      */
#ifdef DEBUG_HKM
      if (dprint) {
        printf(
	    "\tProc %d: Writing global nodes starting at %d to file\n",
	    ProcID, gNodeGlobalMin);
      }
#endif
      if ((fptr = fopen(fname, "ab")) == NULL) {
	(void) fprintf(stderr,"Proc %d: Couldn't open %s: %s\n", ProcID, fname,
                       strerror(errno));
	exit(-1);
      }
      xdrstdio_create(xdrs, fptr, XDR_ENCODE);
      
      /*
      *   Loop over nodes -> This processor will try to process as many
      *   nodes as have contiguous global node numbers, now that it owns
      *   the serial file
      */
      
      while ((gNodeGlobalMin == gNodeLocalNext) && 
             (gNodeGlobalMin < gNum_Nodes)) {
	
	/*
	*   Store the pointer to the node structure corresponding to 
	*   the correct global node
	*/
	
	node = Nodes[oindex];
	
	/*
 	 *   Check to see if this processor actually owns that node
	 */
#ifdef DEBUG_HKM  
	if (! node->Type.Owned) {
	  printf("Error in logic or something\n");
	  exit(-1);
	}
#endif
	/*
	 *    Calculate the block row number of the matrix
	 */
	iblk_row = oindex;
#ifdef DEBUG_HKM
	if (node->Proc_Node_Num != iblk_row) {
	  printf("Screwed up Mesh structure\n");
	  exit(-1);
	}
#endif
	/* 
	*  Value of the beginning row number corresponding to the current
	*  block row. This is also the equation # corresponding to the first
	*  unknown at the node corresponding to the block row.
        *  -> Let's compare against First_Unknown[pg->imtrx] in the node structure
	*     to be sure everything is OK.
	*/
	
	ieqn_row = rpntr[iblk_row];
#ifdef DEBUG_HKM
        if (ieqn_row != node->First_Unknown[pg->imtrx]) {
	  printf("my logic or node_struct is messed up\n");
	  exit(-1);
	}
#endif
	/*
	*  m1 = number of rows in the current row block 
	*/
	
	m1 = rpntr[iblk_row+1] - rpntr[iblk_row];
	
	/*
	*  Cross check the number of unknowns determined above with
	*  the value storred in the var_struct associated with the
	*  node struct
	*/
	
	num_unknowns_node = (int) node->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns;
#ifdef DEBUG_HKM
	if (m1 != num_unknowns_node) {
	  fprintf(stderr,"%s logic error: m1 and ", yo);
	  fprintf(stderr,"num_unknowns_node differ %d %d\n",
	           m1, num_unknowns_node);
	  exit(-1);
	}
#endif
	/*
	*   Calculate the starting index of current row block
	*/
	
	ival = indx[bpntr[iblk_row]];
        a_tmp = a + ival;
	
	/*
	*  Determine the total number of nonzero block columns for the 
	*  current block row and allocate memory for the upcoming sort.
	*/
	
	num_cols = bpntr[iblk_row+1] - bpntr[iblk_row];
	gNodeS = (int *) array_alloc(1, 2*num_cols, sizeof(int));
	iorder = gNodeS + num_cols;
	
	/*
	*  Loop over all the block columns in the current block row. 
	*  Save the global node numbers of the block columns.
	*/
	
	k = 0;
	for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
	  jblk = bindx[j];
	  nodeC = Nodes[jblk];
	  gNodeS[k] = nodeC->Global_Node_Num;
	  iorder[k] = k;
	  k++;
	}
	
	/*
	 *  Sort the global node numbers of the block columns in increasing
	 *  order
	 */
	
	sort2_int_int(num_cols, gNodeS-1, iorder-1);
	
	/*
	 *   Print out Header information for the Current Block Row
	 */
#ifdef DEBUG_HKM
	if (dprint > 2) {
	  fprint_line(stdout, "=", 80);
	  (void) printf("Proc: %d Block Row: %d Local Row Number: %d"
			" Num Rows: %d\n", ProcID, iblk_row, rpntr[iblk_row],
			m1);
	  (void) printf("\t num_block_cols = %d\n", num_cols);
	  fprint_line(stdout, "=", 80);
	}
#endif	
	/*
	*  Dump out Header information for the current block row
	*/
	
	/*       Output the global node number */
        CHECK_XDR(xdr_int(xdrs, &gNodeLocalNext));
	
	/*       Output the processor number */
	CHECK_XDR(xdr_int(xdrs, &ProcID));  
	
	/*       Output the local node number */
	CHECK_XDR(xdr_int(xdrs, &iblk_row));

        /*       Output the coordinates of the current block row
	 *	 node
	 */
        flt_coord[0] = (float) exo->x_coord[iblk_row];
        flt_coord[1] = (float) exo->y_coord[iblk_row];
	if (Num_Dim > 2) {
          flt_coord[2] = (float) exo->z_coord[iblk_row];
	} else {
	  flt_coord[2] = (float) 0.0;
	}
	CHECK_XDR(xdr_vector(xdrs, (char *) flt_coord,  (u_int) 3, 
                             sizeof(float), (xdrproc_t) xdr_float));  
	
	/*       Output the number of unknowns at this node */
        byte_tmp = INT_TO_UCHAR(num_unknowns_node);
	CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));

	/*
	*    Retrieve the initial equation number for the row and store it 
	*    in ieqn.
	*/
	ieqn = ieqn_row;

        /*
         *   Output the first equation number at the current node
         *  -> assume all equation numbers  are contiguous after this one
         */
	CHECK_XDR(xdr_int(xdrs, &ieqn_row));
		
	/*
	 *    Loop over all of the equations at the node, finding the
	 *    variable type, sub_variable type, and checking the node number.
	 */
	for (i = 0; i < num_unknowns_node; i++) {
	  vdT = Index_Solution_Inv(ieqn, &inode, &i_Var_Desc,
		       	           &i_offset, &idof, pg->imtrx);
	  nodeT = Nodes[inode];
#ifdef DEBUG_HKM
	  if (nodeT != node) {
	    fprintf(stderr,"%s Logic error, Proc = %d\n", yo, ProcID);
	    fprintf(stderr,"\t node and nodeT differ %d %d\n",
		    node->Proc_Node_Num, nodeT->Proc_Node_Num);
	    exit (-1);
	  }
#endif
	  /*
	   * Output specific information about the variable
	   */
	  short_tmp = (short int) vdT->Variable_Type;
	  CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	  short_tmp = vdT->MatID;
	  CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	  byte_tmp = INT_TO_UCHAR(vdT->Subvar_Index);
	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
	  byte_tmp = INT_TO_UCHAR(idof);
	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));

	  /*
	   * HKM -> Weights go here. But, will just use a placeholder
	   *        until further research is done
	   */
          flt_tmp = (float) (1.0E-3 * fabs(x[ieqn]) + 1.0E-8);
  	  CHECK_XDR(xdr_float(xdrs, &flt_tmp));
  	  CHECK_XDR(xdr_double(xdrs, x+ieqn));
	  ieqn++;
	}
	
	/*      Output the number of non-zero block columns at this node */
        ushort_tmp = (u_short) num_cols;
        CHECK_XDR(xdr_u_short(xdrs, &ushort_tmp));
	
	/*********************************************************************
         *  Process the block columns in increasing global node number order
	 *********************************************************************/
	
	for (icol = 0; icol < num_cols; icol++) {
	  k = iorder[icol];
	  j = bpntr[iblk_row] + k;
	  jblk = bindx[j];        /* Calculate the local node number of the
				    block column = jblk */
	  /*
	   * Store the location of the current block row-column in a
	   * temporary variable
	   */
	  ival = indx[j];
	  a_tmp = a + indx[j];
	  /*
	  *  Store the pointer to the node structure for the node corresponding
	  *  to the current block column, nodeC.
	  */
	  
	  nodeC = Nodes[jblk]; 
#ifdef DEBUG_HKM
	  if (nodeC->Global_Node_Num != gNodeS[icol]) {
	    (void) fprintf(stderr,"%s Logic error\n", yo);
	    (void) fprintf(stderr,
			   "nodeC->Global_Node_Num and gNodeS[icol] differ: %d %d \n",
			   nodeC->Global_Node_Num, gNodeS[icol]);
	    exit (-1);
	  }
#endif
	  
	  /*
	  *  Calculate the equation number for the first unknown 
	  *  in the current block column
	  */
	  
	  ib1 = cpntr[jblk];
	  
	  /* ending point column index of the current block */
	  
	  ib2 = cpntr[jblk+1];
	  
	  /*
	  *   Calculate the number of columns in the current block column
	  *   and cross check.
	  */
	  
	  n1 = ib2 - ib1;
#ifdef DEBUG_HKM
	  if ((int) nodeC->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns != n1) {
	    fprintf(stderr,"%s Logic error\n", yo);
	    fprintf(stderr,
		    "num unknowns differ: %d %d \n",
		    nodeC->Nodal_Vars_Info[pg->imtrx]->Num_Unknowns, n1);
	    exit (-1);
	  }
#endif
	  /*
	  *     Print out an ascii header information
	  *     -> Mostly for debugging purposes.
	  */
#ifdef DEBUG_HKM 
	  if (dprint > 3) {
	    printf("\t Block Row: %d Block Column: %d"
		   " Num Cols: %d\n", iblk_row, jblk, n1);
	  }
	  if (dprint > 5) {
	    fprint_line(stdout,"-", 80);
	    for (i = 0; i < m1; i++) {
	      for (j = 0; j < n1; j++)
		printf("a[%d]: %e ", ival+j*m1+i, a[ival+j*m1+i]);
	      printf("\n");
	    }
	    fprint_line(stdout,"-", 80);
	  }
#endif
	  /*
	  *  Dump out the header information for the node corresponding to the
          *  block column.
	  */
	  
	  /*        Output global node number of the block column */
  	  CHECK_XDR(xdr_int(xdrs, &gNodeS[icol]));
	  
	  /*        Output the local node number of the block column */
  	  CHECK_XDR(xdr_int(xdrs, &jblk));

          /*       Output the coordinates of the node corresponding to the
                   current column */
          flt_coord[0] = (float) exo->x_coord[jblk];
          flt_coord[1] = (float) exo->y_coord[jblk];
	  if (Num_Dim > 2) {
            flt_coord[2] = (float) exo->z_coord[jblk];
	  } else {
	    flt_coord[2] = (float) 0.0;
	  }
          CHECK_XDR(xdr_vector(xdrs, (char *) flt_coord, (u_int) 3, 
                               sizeof(float), (xdrproc_t) xdr_float));  
	  
	  /*        Output the number of unknowns at the column global node */
          byte_tmp = INT_TO_UCHAR(n1);
  	  CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
	  
	  /*
	  *         Output the identity of first unknown at the current
          *         block column
	  */
          ieqn = nodeC->First_Unknown[pg->imtrx];
  	  CHECK_XDR(xdr_int(xdrs, &ieqn));

	  /*
	  *         Output the identity of the unknowns for the current
	  *         block column
	  */
	  for (i = 0; i < n1; i++) {
            vdT = Index_Solution_Inv(ieqn, &inode, &i_Var_Desc,
		                     &i_offset, &idof, pg->imtrx);
	    nodeT = Nodes[inode];
#ifdef DEBUG_HKM
	    if (nodeT != nodeC) {
	      fprintf(stderr,"%s Logic error, Proc = %d\n", yo, ProcID);
	      fprintf(stderr,"\t nodeC and nodeT differ %d %d\n",
		      nodeC->Proc_Node_Num, nodeT->Proc_Node_Num);
	      exit (-1);
	    }
#endif
	    /*
	     * Output specific information about the variable
	     */
	    short_tmp = (short int) vdT->Variable_Type;
	    CHECK_XDR(xdr_short(xdrs,  &short_tmp));
	    short_tmp = vdT->MatID;
  	    CHECK_XDR(xdr_short(xdrs,  &short_tmp));
            byte_tmp = INT_TO_UCHAR(vdT->Subvar_Index);
  	    CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));
            byte_tmp = INT_TO_UCHAR(idof);
  	    CHECK_XDR(xdr_u_char(xdrs, (unsigned char *) &byte_tmp));

	    
	    ieqn++;
	  }
	  
	  /*
	  *  Dump out the matrix for the current block-row block-column
	  */
	  itmp = m1*n1;
          CHECK_XDR(xdr_vector(xdrs, (char *) a_tmp, (u_int) itmp,
                               sizeof(double), (xdrproc_t) xdr_double));

	} /* END for (icol = 0; icol < num_cols - loop over block cols */
	/*********************************************************************/
	/*
	*   Free allocated memory, allocated within the loop over global
	*   node numbers.
	*/
	
	safer_free((void **) &gNodeS);

	/*
	*  Increment the node on this processor to the node with the
	*  next lowest global node number. End conditions are already
        *  taken care of.
	*/
	
	ii++;
	oindex = ordered_index[ii];
	gNodeLocalNext = ordered_gnodes[ii];

        /*
        *  Increment the global minimum node number
        */

        gNodeGlobalMin++;
	
      } /* End of while loop over processing the local node numbers */

      /*
      *  Close the output file on that processor, and thereby flush
      *  any buffered output
      */
#ifdef DEBUG_HKM
     if (dprint) {
       printf(
 "\tProc = %d, Giving up control: gNodeLocalNext = %d, gNodeGlobalMin = %d\n",
              ProcID, gNodeLocalNext, gNodeGlobalMin);
     }
#endif
      xdr_destroy(xdrs);
      CHECK_FCLOSE(fclose(fptr));
      
    } /* END of -> if (gNodeGlobalMin == gNodeLocalNext) { <- */
    
    /*
    *  Go get the next global minimum node number. And, check to see
    *  whether there still is a processor with more nodes to write out.
    *  If there is, go to the top of the routine and start again.
    */

#ifdef DEBUG_HKM
    if (dprint) {
      if (dprint > 1)
	 (void) printf("\t\tProc = %d, gNodeLocalNext = %d\n", ProcID, 
		       gNodeLocalNext);
      print_sync_end(FALSE);
    }
#endif     
    gNodeGlobalMin = gmin_int(gNodeLocalNext);
  } while(gNodeGlobalMin <gNum_Nodes);
  
  /*
  *  Free memory allocated at the beginning of the routine
  */
  
  safer_free((void **) &ordered_gnodes);
  
  /*
  *  Write an informative message to standard output
  */
  if (ProcID == 0) {
    start_time = ut() - start_time;
    (void) printf("\n\t\tmatrix_dump_vbr: finished, dump time =  %g sec\n", 
		  start_time);
    (void) printf("\t\t"); fprint_line(stdout,"=",64);
  }
} /********************** End matrix_dump_vbr *******************************/
#endif



