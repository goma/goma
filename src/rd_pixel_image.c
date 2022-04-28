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

/*******************************************************************************
 * image_tool.c
 *
 * Tool for importing an image (or a pre-processed image) to calculate relevant
 * field, such as porosity, curvature, etc.  Then interpolates
 * and exports these properties to a pre-made ExodusII mesh.
 *
 * Author:  Scott A Roberts (1514), sarober@sandia.gov
 *
 * Created:  8 March, 2010
 *
 * Modified July 2010,       Ethan Secor        ethan.secor@drake.edu
 *                           Kris Tjiptowidjojo tjiptowi@unm.edu
 * Modified December 2010,   P. R. Schunk prschun@sandia.gov
 *        Bring on board GOMA.
 *
 *******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Include C libraries */
#include "std.h"
/* Include exodus libraries for writing out image */
#include "exodusII.h"
/* Include header files */
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_util.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rd_pixel_image.h"
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_mp.h"
#include "sl_amesos_interface.h"
#include "wr_exo.h"

int first_time_fopen = TRUE;

/*** Begin program *************************************************************/
int rd_image_to_mesh(int N_ext, Exo_DB *exo) {

  /* Input and output file names */
  char outfile[] = "map_pix.exoII";
  char txtfile[] = "efv->name[N_ext]"; // This came in through the input deck

  /* Field variables--These will come in through the input deck as well.  Hardwired here  */
  int num_nod_vars = 1;   // Number of field variables
  char *nod_var_names[1]; // Declare variable names array

  /* File handling variables */
  int CPU_word_size = sizeof(float);
  int IO_word_size = 0;
  FILE *txtid;

  /*Exodus handles */
  float exoversion;
  int exoin;
  float *nodal_var_vals;
  int time_step = 1;
  float time_value = 0.0;

  int e_start, e_end;

  /* variables for processing the image data */
  double minsepar,
      separ; // minimum separation, and separation to match data point to element center

  /* Least square fit variables and arrays */
  double *bf_mat, *f_rhs, *x_fit, *Atranspose_f_rhs;
  int *i_map, *j_map;
  int ipix_blkid = 0;

  /* Bookkeeping arrays */
  int *ElemID_data;

  /* Local variables */
  double nodecoor[MDE][DIM];

  /* Integers */
  int err, i, j, si;
  int elem_loc;
  int ilnode, ignode, ielem, ielem_shape;
  int txt_num_pts = 0;
  int icount;

  /* Array of doubles */
  double *f_data;
  double **xyz_data;
  double **xi_data;
  double **elmctrs;

  /* Scalar doubles */
  double xsum, ysum, zsum;

  double time0;

  time0 = ust();
  /* Scalar floats */

  /* Stop if parallel run */
  if (Num_Proc > 1)
    GOMA_EH(GOMA_ERROR, "pixel mapping is not yet available in parallel. Run serial then use the "
                        "mapped exoII file");

  /* Stop if more than one pixel file read. This will be fixed soon (PRS 8/10/2011 */

  if (N_ext > 0)
    GOMA_EH(GOMA_ERROR, "Sorry, you have to read and convert one pixel file at a time right now.");

  /* Turn matID into blockId --Not necessarly the same */

  int ifound = 0;
  for (i = 0; i < exo->num_elem_blocks; i++) {
    if (efv->ipix_matid[N_ext] == exo->eb_id[i]) {
      ipix_blkid = i;
      ifound = 1;
    }
  }

  if (!ifound)
    GOMA_EH(GOMA_ERROR, "Trouble in rd_pixel_image: cannot find blkid");

  /* Sort out interpolations */
  if ((si = in_list(efv->i[N_ext], 0, Num_Interpolations, Unique_Interpolations)) == -1) {
    GOMA_EH(GOMA_ERROR, "Seems to be a problem finding IntegrationMap interpolation for pixels.");
  }

  /*** Open exodus files and initialize *************************/
  strcpy(txtfile, efv->file_nm[N_ext]);

  /* Replicate current input exodus database to prepare for mapped fields, unless it has already
     been done This is the file we will write out the new fields too for future use as fields
  */

  if (first_time_fopen) {
    one_base(exo, Num_Proc);
    wr_mesh_exo(exo, outfile, 0);
    zero_base(exo);
    first_time_fopen = FALSE;
  }

  /* Now open it again and write to it.  */

  exoin = ex_open(outfile, EX_WRITE, &CPU_word_size, &IO_word_size, &exoversion);
  GOMA_EH(exoin, "ex_open in rd_pixel_image.c");

  /* find element centers- sum node coordinates/numnodes, essentially a center of mass
   * this data is stored in an array, "elmctrs", so it does not need to be recalculated for each
   * data point elmctrs[] holds {xcoord, ycoord, zcoord}
   */

  elmctrs = (double **)malloc(exo->num_elems * sizeof(double *));
  for (ielem = 0; ielem < exo->num_elems; ielem++) {
    elmctrs[ielem] = (double *)malloc(3 * sizeof(double));
  }

  e_start = exo->eb_ptr[ipix_blkid];
  e_end = exo->eb_ptr[ipix_blkid + 1];
  for (ielem = e_start; ielem < e_end; ielem++) {
    load_ei(ielem, exo, 0, pg->imtrx);
    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    for (ilnode = 0; ilnode < exo->eb_num_nodes_per_elem[ipix_blkid]; ilnode++) {
      // ignode = connectivity[ielem][ilnode] - 1;
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      xsum += Coor[0][ignode];
      ysum += Coor[1][ignode];
      if (pd->Num_Dim == 3)
        zsum += Coor[2][ignode];
    }

    elmctrs[ielem][0] = xsum / (double)exo->eb_num_nodes_per_elem[ipix_blkid];
    elmctrs[ielem][1] = ysum / (double)exo->eb_num_nodes_per_elem[ipix_blkid];
    if (pd->Num_Dim == 3)
      elmctrs[ielem][2] = zsum / (double)exo->eb_num_nodes_per_elem[ipix_blkid];
  }

  /*** Find the data point location within the element ***********/

  txtid = fopen(txtfile, "r");
  err = fscanf(txtid, "%i", &txt_num_pts);

  xyz_data = (double **)malloc(txt_num_pts * sizeof(double *));
  for (i = 0; i < txt_num_pts; i++) {
    xyz_data[i] = (double *)malloc(DIM * sizeof(double));
  }

  // remove faulty warning GCC 9
  GOMA_ASSERT_ALWAYS(((size_t)txt_num_pts) < PTRDIFF_MAX);

  f_data = (double *)calloc(txt_num_pts, sizeof(double));

  /* Read the data points */
  for (i = 0; i < txt_num_pts; i++) {
    err = fscanf(txtid, "%le %le %le %le", &(xyz_data[i][0]), &(xyz_data[i][1]), &(xyz_data[i][2]),
                 &(f_data[i]));
  }

  fclose(txtid);

  /* Find element location */

  ElemID_data = (int *)calloc(txt_num_pts, sizeof(int));
  e_start = exo->eb_ptr[ipix_blkid];
  e_end = exo->eb_ptr[ipix_blkid + 1];

  for (i = 0; i < txt_num_pts; i++) {
    minsepar = 1.0e+30; // initialize minimum separation
    elem_loc = -1;      // initialize element ID location
    for (ielem = e_start; ielem < e_end; ielem++) {
      separ = sqrt((xyz_data[i][0] - elmctrs[ielem][0]) * (xyz_data[i][0] - elmctrs[ielem][0]) +
                   (xyz_data[i][1] - elmctrs[ielem][1]) * (xyz_data[i][1] - elmctrs[ielem][1]) +
                   (xyz_data[i][2] - elmctrs[ielem][2]) * (xyz_data[i][2] - elmctrs[ielem][2]));

      if (separ < minsepar) {
        minsepar = separ;
        elem_loc = ielem + 1;
      }
    }

    ElemID_data[i] = elem_loc;
  }

  /* Find local coordinates */
  xi_data = (double **)malloc(txt_num_pts * sizeof(double *));
  for (i = 0; i < txt_num_pts; i++) {
    xi_data[i] = (double *)malloc(DIM * sizeof(double));
  }

  /* Here we are going to pull a fast one.  For deforming mesh problems, you must have
   * the elemenet_stiffness_pointers set up to point to the displacements in the solution vector.
   * Trouble is, load_elem_dof_ptrs has not been called yet, and so when beer_belly is called by
   * find_xi below, it bombs.   DeformingMesh is true but *esp is null.   So, we will just shut
   * off deforming mesh for purposes of the map, and then turn it back on if it was already on
   */
  int if_deform = 0;

  for (i = 0; i < txt_num_pts; i++) {

    /* Hold on, I may need load_elem_dofptr here.  and a call to bf_mp_init(pd) as in mm_fill.c */
    ielem = ElemID_data[i] - 1;
    ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);
    load_ei(ielem, exo, 0, pg->imtrx);

    if (ei[pg->imtrx]->deforming_mesh) {
      if_deform = 1;
      ei[pg->imtrx]->deforming_mesh = FALSE;
    }

    for (ilnode = 0; ilnode < exo->eb_num_nodes_per_elem[ipix_blkid]; ilnode++) {
      // ignode = connectivity[ielem][ilnode] - 1;
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      for (j = 0; j < pd->Num_Dim; j++) {
        nodecoor[ilnode][j] = Coor[j][ignode];
      }
    }

    find_xi(ielem, xyz_data[i], xi_data[i], ei[pg->imtrx]->ielem_type, nodecoor, si, N_ext);
    /* Now undo the mess if needed  */
    if (if_deform)
      ei[pg->imtrx]->deforming_mesh = TRUE;
  }

  /*** Assemble matrix for least square fit ***********/
  /* Here is where we may just load one row at a time, viz. allocate bf_mat to just exo->num_nodes
     vector and then call buld_epetra_matrix[row_i}. We are also going to want to skip zero columns
     here
  */

  f_rhs = (double *)malloc(txt_num_pts * sizeof(double));
  x_fit = (double *)calloc(exo->eb_num_nodes_per_elem[ipix_blkid] * exo->eb_num_elems[ipix_blkid],
                           sizeof(double));
  Atranspose_f_rhs = (double *)malloc(exo->eb_num_nodes_per_elem[ipix_blkid] *
                                      exo->eb_num_elems[ipix_blkid] * sizeof(double));
  bf_mat = (double *)malloc(txt_num_pts * exo->eb_num_nodes_per_elem[ipix_blkid] * sizeof(double));
  i_map = (int *)malloc(txt_num_pts * exo->eb_num_nodes_per_elem[ipix_blkid] * sizeof(int));
  j_map = (int *)malloc(txt_num_pts * exo->eb_num_nodes_per_elem[ipix_blkid] * sizeof(int));

  for (i = 0; i < txt_num_pts; i++) {
    f_rhs[i] = f_data[i];
  }

  for (i = 0; i < txt_num_pts * exo->eb_num_nodes_per_elem[ipix_blkid]; i++) {
    bf_mat[i] = 0.0;
    i_map[i] = j_map[i] = -1;
  }

  /* Three steps left */
  /* 1) compute bf_mat and i_map, j_map  (Both arrays are
        0,num_txt_points*num_nodes long
        bf_mat has the matrix values, and i_map is the corresponding row and j_map the column.
        Note that we have compressed out all the zeros this way.
     2) transport bf_mat, i_map, j_map to a C++ routine trilinos_solve_ls
     3) trilinos_solve_ls builds Epetra CSR matrix. Directs normal equation solve.
  */

  /* Step 1 */
  icount = 0;
  for (i = 0; i < txt_num_pts; i++) {
    ielem = ElemID_data[i] - 1;
    err = load_ei(ielem, exo, 0, pg->imtrx);
    GOMA_EH(err, "load_ei");

    ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);
    ielem_shape = type2shape(ei[pg->imtrx]->ielem_type);

    for (ilnode = 0; ilnode < exo->eb_num_nodes_per_elem[ipix_blkid]; ilnode++) {
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      i_map[icount] = i;
      j_map[icount] = ignode;
      // Note hardwired to be Q1 now. But this should come from input deck
      bf_mat[icount] = newshape(xi_data[i], ei[pg->imtrx]->ielem_type, PSI, ilnode, ielem_shape,
                                efv->i[N_ext], ilnode);
      icount++;
    }
  }

  /*** Step 2, 3: transport bf_mat, i_map, and j_map to trilinos and  Solve the least squares system
   * ************/
  int num_blk_nodes = exo->eb_num_nodes_per_elem[ipix_blkid] * exo->eb_num_elems[ipix_blkid];

  trilinos_solve_ls(bf_mat, j_map, f_rhs, x_fit, Atranspose_f_rhs, txt_num_pts,
                    exo->eb_num_nodes_per_elem[ipix_blkid], num_blk_nodes, 1);

  /*** Setup variables in exodus file ******************************************/

  // Number of variables and names
  nod_var_names[0] = efv->name[N_ext];
  err = ex_put_variable_param(exoin, EX_NODAL, num_nod_vars);
  err = ex_put_variable_names(exoin, EX_NODAL, num_nod_vars, nod_var_names);

  // Set time
  err = ex_put_time(exoin, time_step, &time_value);

  nodal_var_vals = (float *)calloc(exo->num_nodes, sizeof(float));
  for (i = 0; i < exo->num_nodes; i++)
    nodal_var_vals[i] = 0.;

  // Loop over all of the nodes in the mesh block of interest
  e_start = exo->eb_ptr[ipix_blkid];
  e_end = exo->eb_ptr[ipix_blkid + 1];
  for (ielem = e_start; ielem < e_end; ielem++) {
    load_ei(ielem, exo, 0, pg->imtrx);
    for (ilnode = 0; ilnode < exo->eb_num_nodes_per_elem[ipix_blkid]; ilnode++) {
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      // printf(" node = %d, val = %f \n", i+1, x_fit[i]);
      nodal_var_vals[ignode] = (float)x_fit[ignode];
    } // End of loop over all of the nodes (i)
  }

  err = ex_put_var(exoin, time_step, EX_NODAL, 1, 1, exo->num_nodes, nodal_var_vals);
  /* Now allocate memory for external field variable array for use and populate */
  /* You need to do this even though you have writtenit to a file */

  efv->ext_fld_ndl_val[N_ext] = alloc_dbl_1(exo->num_nodes, 0.0);
  printf("rd_pix_image: Allocated field %d for %s at %p\n", N_ext, efv->name[N_ext],
         (void *)efv->ext_fld_ndl_val[N_ext]);
  for (i = 0; i < exo->num_nodes; i++) {
    efv->ext_fld_ndl_val[N_ext][i] = nodal_var_vals[i];
  }

  printf("Time for pixel-to-mesh processing: %g seconds\n", ust() - time0);

  /* Calculate error*/
  /* This uses the same function as in rd_pixel_image2.c, so we need to do a few
     'conversions' to accomodate the different pixel data structure the function expects
   */
  double resolution[3];
  /*   double pixorigin[3]; */
  int pixsize[3];
  double ***pixdata;
  double minx, miny, maxx, maxy;
  double firstx = 0;
  int setres;
  double dx = 0, dy = 0;
  int ii, jj;
  // double rmsd_error;
  double pmax, pmin;
  minx = miny = 1e20;
  maxx = maxy = -1e20;
  setres = 0;
  resolution[0] = 0;
  resolution[1] = 0;
  for (i = 0; i < txt_num_pts; i++) {
    if (i == 0)
      firstx = xyz_data[i][0];
    if (xyz_data[i][0] < minx)
      minx = xyz_data[i][0];
    if (xyz_data[i][0] > maxx)
      maxx = xyz_data[i][0];
    if (xyz_data[i][1] < miny)
      miny = xyz_data[i][1];
    if (xyz_data[i][1] > maxy)
      maxy = xyz_data[i][1];
    if (xyz_data[i][0] != firstx && !setres) {
      resolution[0] = resolution[1] = dx = dy = xyz_data[i][0] - firstx;
      setres = 1;
    }
  }

  /*   pixorigin[0] = minx;
       pixorigin[1] = miny;
       pixorigin[2] = 0;
  */
  pixsize[0] = (int)((maxx - minx) / resolution[0]) + 1;
  pixsize[1] = (int)((maxy - miny) / resolution[1]) + 1;
  pixsize[2] = 1;
  pixdata = (double ***)malloc(pixsize[0] * pixsize[1] * pixsize[2] * sizeof(double **));
  for (i = 0; i < pixsize[0]; i++) {
    pixdata[i] = (double **)malloc(pixsize[1] * pixsize[2] * sizeof(double *));
    for (j = 0; j < pixsize[1]; j++) {
      pixdata[i][j] = (double *)malloc(pixsize[2] * sizeof(double));
    }
  }
  pmax = -1e20;
  pmin = 1e20;
  for (i = 0; i < txt_num_pts; i++) {
    ii = (int)((xyz_data[i][0] - minx) / dx);
    jj = (int)((xyz_data[i][1] - miny) / dy);
    pixdata[ii][jj][0] = f_data[i];
    if (f_data[i] > pmax)
      pmax = f_data[i];
    if (f_data[i] < pmin)
      pmin = f_data[i];
  }

  // printf("Calculating error...\n");
  // All that just to run calc_error..
  // rmsd_error = calc_error(pixdata, pixsize, resolution, pixorigin, nodal_var_vals, exo,
  // ipix_blkid, si, N_ext); printf("Scaled RMSD error was %g %%\n",rmsd_error / (pmax - pmin));

  /* Now you can clean things up */
  err = ex_close(exoin);

  safe_free(pixdata);
  safe_free(elmctrs);
  safe_free(xyz_data);
  safe_free(f_data);
  safe_free(ElemID_data);
  safe_free(xi_data);
  safe_free(f_rhs);
  safe_free(x_fit);
  safe_free(Atranspose_f_rhs);
  safe_free(bf_mat);
  safe_free(i_map);
  safe_free(j_map);
  safe_free(nodal_var_vals);

  return (0);
}

/******************************************************************************/
/* Note that this routine is a simpler version of two others in Goma which do things that we don't
 * need to here.  Specifically, cf. with get_element_xi_newton and invert_isoparametric_map()
 */

extern int find_xi(int elem_id,               /*known element id number*/
                   const double x[DIM],       /*x,y,z coordinates of data point*/
                   double xi[DIM],            /*local coordinate output*/
                   int elem_type,             /*element type*/
                   double nodecoor[MDE][DIM], /*global coordinates of local nodes  nodecoor[local
                                                 node number][x y or z]*/
                   int si,
                   int N_ext) /*Interpolation index */
{

  /* local variables*/
  int i, j, error;
  int ielem_shape;
  int nodenum;
  int converged, inewton; /*flag for convergence and counter for iterations*/

  double R[DIM];              // residual vector
  double update[DIM] = {0.0}; // update for xi[]

  double norm; // convergence norm

  double phi[MDE]; // basis functions, one for each local node

  double x_intp[DIM]; // sum over local nodes of global coordinates

  converged = 0;
  inewton = 0;

  nodenum = elem_info(NNODES, elem_type);
  ielem_shape = type2shape(elem_type);

  xi[0] = xi[1] = xi[2] = 0.;
  while ((!converged) && (inewton < 50)) {
    error = load_basis_functions(xi, bfd);
    GOMA_EH(error, "load_basis_functions");

    error = beer_belly();
    GOMA_EH(error, "beer_belly");

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      x_intp[i] = 0.0;
      for (j = 0; j < nodenum; j++) {
        phi[j] = newshape(xi, elem_type, PSI, j, ielem_shape, efv->i[N_ext], j);
        x_intp[i] += nodecoor[j][i] * phi[j];
      }
      // x[i] is global xyz coordinate of the point, x_intp[i] is sum over local node basis
      // functions of global coordinates

      R[i] = x[i] - x_intp[i];
      // the residual
    }

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      update[i] = 0.0;
    }
    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        update[i] += bfd[si]->B[j][i] * R[j];
      }
    }

    norm = 0.0;
    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      norm += update[i] * update[i];
    }
    converged = (norm < 1.e-10);

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      xi[i] += update[i];
    }

    inewton++;
  }
  return (converged);
}

/*** end program *************************************************************/
