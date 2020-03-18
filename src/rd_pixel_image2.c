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

/******************************************************************************* 
 * rd_pixel_image2.c
 *
 * Tool for mapping a pixel or voxel file (2D or 3D image, respectively) 
 * to material property (e.g. thermal/electrical conductivity,
 * porosity, shell height function, etc.). Produces a nodal field similar 
 * to an 'External Field' that can be used as input for material properties. 
 * The field is saved as an ExodusII file called 'map_pix_fast.exoII'.
 *
 *
 * Authors: Dan S Bolintineanu (1516) dsbolin@sandia.gov, P.R. Schunk (1516) prschun@sandia.gov
 *
 * 
 * Created:  Sept. 16, 2013
 *
 *
 *******************************************************************************/

/* Include C libraries */
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>

#define GOMA_MM_FILL_C
#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_particles.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_update_parameter.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc_rotate.h"
#include "bc_special.h"
#include "bc_surfacedomain.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "loca_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rd_pixel_image.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_amesos_interface.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "std.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"


int first_time_fopen2 = TRUE;

static double bi_interp
(double xg, double yg, 
	double resolution[],
	double ***pixdata,
	int pixsize[]);

static double tri_interp
(double xg, double yg, double zg, 
	double resolution[],
	double ***pixdata,
	int pixsize[]);

/*** Begin program *************************************************************/
int
rd_image_to_mesh2(int N_ext, Exo_DB *exo)  
{

  /* Input and output file names */
  char exooutfilename[256];   
  char voxfilename[256];

  /* Field variables.  */
  int num_nod_vars;                  // Number of field variables (no duplicates)
  int dupflag;                       // Flag for checking duplicates
  char ** nod_var_names;             // Variable names array, no duplicates
  char *nod_var_names_temp[efv->Num_external_pixel_field]; // Variable names array, including possible duplicates
 
  /* File handling variables */
  int CPU_word_size = sizeof(float);
  int IO_word_size = 0;
  FILE * pixfid;

  /* Pixel file descriptors */  
  int ipix_blkid = 0;           // Block ID to be filled
  int pixdim;                   // Dimensionality (2 or 3)
  int pixsize[3];               // Size in each dimension
  double ***pixdata;            // All data points
  double pmax, pmin;            // Max/min pixel/voxel value
  int numpix;                   // Number of data points
  double resx, resy, resz;      // Size of pixels/voxels in x,y,z directions
  double resolution[3];         // Array form of above (easier to pass to function)
  int pixinterp = 0;            // Flag indicating pixel interpolation (0 - none; 1 - bi/trilinear)
  double x0, y0, z0;            // Coordinates of origin for pixel/voxel file
  /*  double pixorigin[3];           // Array form of above */

  float exoversion;
  int exoout;
  float *nodal_var_vals;
  float *nodal_var_vals_err;
  int time_step = 1;
  double time_value = 0.0;

  int e_start, e_end;  
  int curr_eb_nodes;            // Total number of nodes in elem block 
  int *my_ignodes;              // List of global node ids for elem block

  /*Lumped mass variables*/
  double * lumped_mass;
  double * bvec;
  double *save_gp;
  int gptot, gp_global;

  double xi[3];

  int err, i, j, k, si;
  int igp, ngp;
  int ilnode, jlnode, ignode, ielem;  
  int ix, iy, iz;

  double xg, yg, zg;
  double phi_i, phi_j;
  double val;
  double wt, det_J;

  int nmiss;
  //  double rmsd_error;
  double time0;

  /*  int first_time_field = TRUE; // Flag to check whether this is the first time mapping to the specified field */

  int my_N_ext; // Index of external field being mapped to. May not be what's passed to the function.

  time0 = ust();
  my_N_ext = N_ext;

  /* Stop if parallel run */
  if(Num_Proc > 1) EH(-1, "pixel mapping is not yet available in parallel. Run serial then use the mapped exoII file(s)");

  /* Turn matID into blockId --Not necessarly the same */
  int ifound = 0;
  for (i=0; i<exo->num_elem_blocks; i++) {
    if(efv->ipix_matid[N_ext] == exo->eb_id[i]){
      ipix_blkid = i;
      ifound = 1;
    }
  }

  if(!ifound) EH(-1,"Trouble in rd_pixel_image: cannot find blkid");

  curr_eb_nodes = exo->eb_num_elems[ipix_blkid]*exo->eb_num_nodes_per_elem[ipix_blkid];
  my_ignodes = (int *) malloc(curr_eb_nodes * sizeof(int));

  /* Sort out interpolations */
  if (( si = in_list(efv->i[N_ext], 0, Num_Interpolations, 
		     Unique_Interpolations)) == -1) {
    EH(-1,"Seems to be a problem finding IntegrationMap interpolation for pixels.");
  }


  /***************************/
  /* Open exodus output file */
  /***************************/
  sprintf(exooutfilename,"%s.exoII","map_pix_fast");

  /* Replicate current input exodus database to prepare for mapped fields, unless it has already been done
     This is the file we will write out the new fields to for future use as fields */

  if (first_time_fopen2){
    one_base(exo);
    wr_mesh_exo(exo, exooutfilename, 0); 
    zero_base(exo);
    
    exoout = ex_open(exooutfilename, EX_WRITE, &CPU_word_size, &IO_word_size, &exoversion);
    EH(exoout, "ex_open in rd_pixel_image2.c");

    num_nod_vars = 0;
    for (i = 0; i < efv->Num_external_field; i++){
      if (efv->ipix[i] == 2){
	dupflag = 0;
	for (j = 0; j < num_nod_vars; j++){
	  if (strcmp(nod_var_names_temp[j], efv->name[i]) == 0){
	    dupflag = 1;
	  }
	}
	if (!dupflag){
	  nod_var_names_temp[num_nod_vars++] = efv->name[i];
	}		     
      }
    }
    nod_var_names = (char **) malloc(num_nod_vars * sizeof(char*));
    for (i = 0; i < num_nod_vars; i++){
      nod_var_names[i] = (char *)malloc(MAX_STR_LENGTH * sizeof(char));
      nod_var_names[i] = nod_var_names_temp[i];
    }
    
    err = ex_put_variable_param(exoout, EX_NODAL, num_nod_vars);
    err = ex_put_variable_names(exoout, EX_NODAL, num_nod_vars, nod_var_names);
    err = ex_put_time(exoout, time_step, &time_value);
    safe_free(nod_var_names);
    first_time_fopen2 = FALSE;
  }
  else{
    exoout = ex_open(exooutfilename, EX_WRITE, &CPU_word_size, &IO_word_size, &exoversion);
    EH(exoout, "ex_open in rd_pixel_image2.c");
    //Loop through all fields, check for same variable being mapped from another field
    for (i = N_ext-1; i >= 0; i--){
      if (!strcmp(efv->name[i],efv->name[N_ext])){
	my_N_ext = i;
	/*	first_time_field = FALSE; */
	break;
      }
    }
  }

  /****************************************/
  /* Open and read voxel file             */
  /****************************************/
  sprintf(voxfilename,"%s",efv->file_nm[N_ext]);
  pixfid = fopen(voxfilename, "r"); 
  if (pixfid == NULL) {
    EH(-1,"Could not open voxel file!");
  }
  err = fscanf(pixfid,"%d", &pixdim);

  printf("Reading pixel file %s of dimensionality %d into field %s\n",voxfilename,pixdim,efv->name[my_N_ext]);

  pixsize[2] = 1;
  resz = 1;
  z0 = 0;
  if (pixdim == 2) {
    err = fscanf(pixfid,"%d %d",&(pixsize[0]),&(pixsize[1]));
    if (err != 2) {
      EH(-1, "Error reading pixel file expected two ints");
      return -1;
    }
    err = fscanf(pixfid,"%lf %lf",&resx, &resy);
    if (err != 2) {
      EH(-1, "Error reading pixel file expected two floats");
      return -1;
    }
    err = fscanf(pixfid,"%lf %lf",&x0, &y0);
    if (err != 2) {
      EH(-1, "Error reading pixel file expected two floats");
      return -1;
    }
  }
  else if (pixdim == 3){
    err = fscanf(pixfid,"%d %d %d",&(pixsize[0]),&(pixsize[1]),&(pixsize[2]));
    if (err != 3) {
      EH(-1, "Error reading pixel file expected three ints");
      return -1;
    }
    err = fscanf(pixfid,"%lf %lf %lf",&resx, &resy, &resz);
    if (err != 3) {
      EH(-1, "Error reading pixel file expected three floats");
      return -1;
    }
    err = fscanf(pixfid,"%lf %lf %lf",&x0, &y0, &z0);
    if (err != 3) {
      EH(-1, "Error reading pixel file expected three floats");
      return -1;
    }
  }
  else {
    EH(-1,"Problem reading pixel/voxel input file; first entry should be dimensionality (2 or 3)");
    return -1;
  }

  resolution[0] = resx;
  resolution[1] = resy;
  resolution[2] = resz;

  /*  pixorigin[0] = x0;
      pixorigin[1] = y0;
      pixorigin[2] = z0;
  */

  if (pixdim != pd->Num_Dim) {
    printf("WARNING: Pixel file dimensionality is %d, problem dim. is %d\n",pixdim,pd->Num_Dim);
  }

  numpix = pixsize[0]*pixsize[1]*pixsize[2];
  
  pixdata = (double ***) malloc(numpix * sizeof(double **));
  for (i = 0; i < pixsize[0]; i++)
    {
      pixdata[i] = (double **) malloc(pixsize[1]*pixsize[2]*sizeof(double*));
      for (j = 0; j < pixsize[1]; j++)
	{
	  pixdata[i][j] = (double *) malloc(pixsize[2]*sizeof(double));
	}
    }

  /* Read the data points */
  for (i = 0; i < pixsize[0]; i++)
    for (j = 0; j < pixsize[1]; j++)
      for (k = 0; k < pixsize[2]; k++)
	{
	  err = fscanf( pixfid, "%lf ", &(pixdata[i][j][k]));
	  EH(err,"Problem reading file!\n");
	}

  fclose(pixfid);


  /****************************************/
  /* Some initialization things..         */
  /****************************************/


  lumped_mass = (double *) malloc(exo->num_nodes * sizeof(double));
  bvec = (double *) malloc(exo->num_nodes * sizeof(double));

  for (i = 0; i < exo->num_nodes; i++)
    {
      lumped_mass[i] = bvec[i] = 0;
    }

  e_start = exo->eb_ptr[ipix_blkid]; e_end = exo->eb_ptr[ipix_blkid+1];

  gptot = gp_global = 0;
  j = 0;
  //Loop over elements; count GQ points, store list of global node ids
  for (ielem = e_start; ielem < e_end; ielem++){
    load_ei(ielem, exo, 0, pg->imtrx); 
    ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);
    gptot += elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
    for (ilnode = 0; ilnode < ei[pg->imtrx]->num_local_nodes; ilnode++) {                                           
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      my_ignodes[j++] = ignode;
    }
  }

  save_gp = (double *) malloc(gptot * sizeof(double));

  for (i = 0; i < gptot; i++) {
    save_gp[i] = efv->empty_value[N_ext];
  }

  /*************************************************************/
  /* Loop over elements to assign pixel values to Gauss points */
  /*************************************************************/ 

  nmiss = 0;
  //Loop over elements to assign pixel values to Gauss points
  for (ielem = e_start; ielem < e_end; ielem++) 
    {
      load_ei(ielem, exo, 0, pg->imtrx);
      ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);

      err = bf_mp_init(pd);
      EH(err, "bf_mp_init");
      
      ngp = elem_info(NQUAD, ei[pg->imtrx]->ielem_type); //Number of Gauss points (GP) in element
      
      for (igp = 0; igp < ngp; igp++) //Loop over GP's
	{ 	  
	  err = load_basis_functions(xi, bfd);
	  EH( err, "problem from load_basis_functions");
	  
	  err = beer_belly();
	  EH( err, "beer_belly");
	  
	  find_stu(igp, ei[pg->imtrx]->ielem_type, &(xi[0]), &(xi[1]), &(xi[2]));
	  
	  xg = yg = zg = 0.0;
	  for (ilnode = 0; ilnode < ei[pg->imtrx]->num_local_nodes; ilnode++) {
	    ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode]; 
	    phi_i = newshape(xi, ei[pg->imtrx]->ielem_type, PSI, ilnode, ei[pg->imtrx]->ielem_shape, efv->i[N_ext], ilnode);
	    xg += Coor[0][ignode]*phi_i;
	    yg += Coor[1][ignode]*phi_i;
	    if (pixdim == 3) zg += Coor[2][ignode]*phi_i;
	  }
	  //xg, yg, zg are now the global coordinates of the Gauss point
	  //Get the pixel value at this point
	  if (pixinterp){	    
	    if (pixsize[2] == 1) { //bilinear interpolation
	      val = bi_interp(xg,yg,resolution,pixdata,pixsize);
	    }
	    else { //trilinear interpolation
	      val = tri_interp(xg,yg,zg,resolution,pixdata,pixsize);
	    }
	  }
	  else {
	    ix = (int) ((xg - x0) / resx);
	    iy = (int) ((yg - y0) / resy);
	    iz = (int) ((zg - z0) / resz);
	    if (ix < 0 || ix >= pixsize[0] ||
		iy < 0 || iy >= pixsize[1] ||
		iz < 0 || iz >= pixsize[2]){
	      nmiss++;
	      val = efv->empty_value[N_ext];
	    }
	    else {
	      val = pixdata[ix][iy][iz]; 
	    }
	  }
	  save_gp[gp_global++] = val;
	}
    }

  if (nmiss){
    printf("WARNING: Some portions of the mesh/element block being mapped to are located outside the range of the voxel field (%d of %d Gauss points)\n",nmiss,gp_global);
  }
 
  gp_global = 0;
  //printf("Gauss point values assigned\n");

  /***************************************************************/
  /* Now loop over elements for lumped mass least square method  */
  /***************************************************************/ 

  for (ielem = e_start; ielem < e_end; ielem++) //Loop over elements
    {
      load_ei(ielem, exo, 0, pg->imtrx);
      ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);

      err = bf_mp_init(pd);
      EH(err, "bf_mp_init");

      for (igp = 0; igp < elem_info(NQUAD, ei[pg->imtrx]->ielem_type); igp++)
	{
	  err = load_basis_functions( xi, bfd);
	  EH( err, "problem from load_basis_functions");

	  err = beer_belly();
	  EH( err, "beer_belly");
	  
	  det_J = bfex[N_ext]->detJ;
	  wt = Gq_weight(igp, ei[pg->imtrx]->ielem_type);
       	  
	  val = save_gp[gp_global++];

	  for (ilnode = 0; ilnode < ei[pg->imtrx]->num_local_nodes; ilnode++) 
	    {
	      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
	      phi_i = bfex[N_ext]->phi[ilnode];
	      bvec[ignode] += phi_i * val * wt * det_J;
	      for (jlnode = 0; jlnode < ei[pg->imtrx]->num_local_nodes; jlnode++)
		{
		  phi_j = bfex[N_ext]->phi[jlnode];
		  lumped_mass[ignode] += phi_i*phi_j*wt*det_J;
		}
	    }
	}
    }

  
  /*************************************/
  /* Output to map_pix exoII file      */
  /* and fill in ext_fld_ndl_val array */
  /*************************************/

  nodal_var_vals_err = (float *) calloc(curr_eb_nodes, sizeof(float));
  nodal_var_vals = (float *) calloc(exo->num_nodes, sizeof(float));
  efv->ext_fld_ndl_val[N_ext] = alloc_dbl_1(exo->num_nodes, 0.0); 

  for (i = 0; i < curr_eb_nodes; i++){ 
    nodal_var_vals_err[my_ignodes[i]] = 
      efv->ext_fld_ndl_val[my_N_ext][my_ignodes[i]] = 
      bvec[my_ignodes[i]]/lumped_mass[my_ignodes[i]];
  }
  
  for (i = 0; i < exo->num_nodes; i++){
    nodal_var_vals[i] = efv->ext_fld_ndl_val[my_N_ext][i];
  }

  k = 0;
  for (i = 0; i < efv->Num_external_field; i++){
    if (efv->ipix[i] == 2){
      k++;
    }
    if (!(strcmp(efv->name[i], efv->name[my_N_ext]))){
      break;
    }
  }

  err = ex_put_var(exoout, time_step, EX_NODAL, k, 1, exo->num_nodes, nodal_var_vals);



  /*****************************/
  /* Timing and error estimate */
  /*****************************/

  printf("Time for pixel-to-mesh processing: %6.4lf seconds\n",ust()-time0);

  pmin = 1e20;
  pmax = -1e20;
  for (i = 0; i < pixsize[0]; i++)
    for (j = 0; j < pixsize[1]; j++)
      for (k = 0; k < pixsize[2]; k++){
	if (pixdata[i][j][k] < pmin) pmin = pixdata[i][j][k];
	if (pixdata[i][j][k] > pmax) pmax = pixdata[i][j][k];	
      }
  
  /*  rmsd_error = calc_error(pixdata, pixsize, resolution, pixorigin, nodal_var_vals_err, exo, ipix_blkid, si, my_N_ext);
  printf("Scaled RMSD error was %g %%\n",rmsd_error/(pmax-pmin)*100.0);
  if (pixsize[2] > 1) {
    printf("Warning: For 3D problems with unstructured or tet meshes, this error value may be inaccurate (tends to be overestimated)\n");    
  }
  */
  
  /***********************/
  /* Clean everything up */
  /***********************/

  safe_free(nodal_var_vals_err);
  safe_free(pixdata);
  safe_free(lumped_mass);
  safe_free(bvec);
  safe_free(save_gp);
  safe_free(my_ignodes);

  return(0);
}


extern double calc_error(double ***pixdata, 
			 int pixsize[], 
			 double resolution[],
			 double pixorigin[],
			 float *nodal_var_vals, 
			 Exo_DB *exo, 
			 int ipix_blkid,
			 int si,
			 int N_ext)
{
  /**************************************************************/
  /* Routine to calculate the error as follows: */

  /* - Loop through all elements. */
  /*  - Loop through all pixels to find those that belong to each element. */
  /*    This is not entirely trivial - the current implementation */
  /*    uses a simple bounding box check as a first pass, followed by a ray-casting algorithm. There */
  /*    are probably easier ways, but this works without fail for any polygonal 2D element.  */
  /*  - For each such pixel, calculate the interpolated value at that pixel's center, compare */
  /*    to the actual value, and add the square difference to the total error. */

  /* Note: The current implementation of this error calculation is rigorous for 2D only  */
  /* (mostly due to issues with figuring out which element a voxel belongs to in 3D).  */
  /* In the original pixel-to-mesh tool, something similar is done at some point.  */
  /* It involves an assumption that each pixel belongs to the element with  */
  /* the closest center, which is not true in general. For 3D, we just use a bounding box check here,  */
  /* which is arguably worse than the assumption in the original pixel-to-mesh code (it wouldn't be */
  /* a big deal to also do the 'closest-center' check, but it doesn't really 'mesh' */
  /* with the loop structure here). */

  /* Both the bounding box and closest-center checks are problematic for error  */
  /* calculation purpposes, particularly for tet. meshes or  */
  /* unstructured meshes in general. Then again, why would you be using anything  */
  /* other than a structured quad/hexahedral mesh for voxel-to-mesh problems?  */

  /* A warning message is printed to this effect for 3D situations. */
  
  double tot_error = 0;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  int local_elem_node_id[MAX_NODES_PER_SIDE];
  int nodes_per_side;
  int ncoll;
  double px, py, pz, rx, ry;
  double n1x, n1y, n2x, n2y;
  double ux, uy;
  double n1p_x, n1p_y;
  double r_cross_u, n1p_cross_u, n1p_cross_r;

  double sparam, tparam;
  double nodecoor[MDE][DIM];
  double pixcoords[3];
  double xi_data[3];

  double phi_i;
  double pix_interp;
  int i, j, k, iside, ilnode, ignode, ielem;
  int e_start, e_end;

  int max_error_points = 2000; 
  //Probably should have some better heuristic for this that takes into account
  //number of voxels..


  printf("Calculating error..\n");
  //For ray-casting:
  rx = 0.13215513585813;
  ry = sqrt(1-rx*rx);

  int numpoints = 0;

  e_start = exo->eb_ptr[ipix_blkid]; e_end = exo->eb_ptr[ipix_blkid+1]; 

  for (ielem = e_start; ielem < e_end; ielem++) {
    if (e_end - e_start > max_error_points){
      ielem += (e_end - e_start)/max_error_points;
    }
    load_ei(ielem, exo, 0, pg->imtrx);
    ei[pg->imtrx]->ielem_type = Elem_Type(exo, ielem);
    //Get bounding box  
    xmin = ymin = zmin = 1e20;
    xmax = ymax = zmax = -1e20;
    	  
    for (ilnode = 0; ilnode < ei[pg->imtrx]->num_local_nodes; ilnode++) {
      ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
      if (Coor[0][ignode] < xmin) xmin = Coor[0][ignode];
      if (Coor[0][ignode] > xmax) xmax = Coor[0][ignode];
      if (Coor[1][ignode] < ymin) ymin = Coor[1][ignode];
      if (Coor[1][ignode] > ymax) ymax = Coor[1][ignode];
      nodecoor[ilnode][0] = Coor[0][ignode];	     
      nodecoor[ilnode][1] = Coor[1][ignode];
      nodecoor[ilnode][2] = Coor[2][ignode];
      if (pixsize[2] > 1){
	 if (Coor[2][ignode] < zmin) zmin = Coor[2][ignode];
	 if (Coor[2][ignode] > zmax) zmax = Coor[2][ignode];
      }    
    }

    for (i = 0; i < pixsize[0]; i++){
      px = pixorigin[0]+(i+0.5)*resolution[0];
      for (j = 0; j < pixsize[1]; j++){
	py = pixorigin[1]+(j+0.5)*resolution[1];
	for (k = 0; k < pixsize[2]; k++){
	  pz = pixorigin[2]+(k+0.5)*resolution[2];
	  if (pixsize[2] == 1) {
	    pz = pixcoords[2] = pixorigin[2];
	  }

	  // Bounding box test
	  if (px < xmin || px > xmax || py < ymin || py > ymax) continue;
	  if (pixsize[2] > 1){
	    if (pz < zmin || pz > zmax) {
	      continue;
	    }	  
	    else{
	      pixcoords[0] = px;
	      pixcoords[1] = py;
	      pixcoords[2] = pz;
	    }
	  }	  
	  else{
	  // Ray collision test:
	    ncoll = 0;
	    for (iside = 0; iside < ei[pg->imtrx]->num_sides; iside++){	 
	      get_side_info(ei[pg->imtrx]->ielem_type, iside+1, &nodes_per_side, local_elem_node_id);
	      
	      n1x = Coor[0][Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + local_elem_node_id[0]]];
	      n1y = Coor[1][Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + local_elem_node_id[0]]];
	      
	      n2x = Coor[0][Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + local_elem_node_id[1]]];
	      n2y = Coor[1][Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + local_elem_node_id[1]]];

	      ux = n2x - n1x;
	      uy = n2y - n1y;
	      
	      //Check this.. it may not be right yet..
	      ux*= 10*resolution[0];
	      uy*= 10*resolution[1];

	      r_cross_u = (rx - uy)*(ry - uy);	    
	      if (r_cross_u == 0) continue; //Unlikely case that ray is perfectly parallel to an element edge
	      n1p_x = n1x - px;
	      n1p_y = n1y - py;
	      
	      n1p_cross_u = (n1p_x - uy)*(n1p_y - ux);
	      sparam = n1p_cross_u / r_cross_u;
	      if (sparam < 0) continue;
	      
	      n1p_cross_r = (n1p_x - ry)*(n1p_y - rx);
	      tparam = n1p_cross_r/r_cross_u;
	      if (tparam < 0 || tparam > 1) continue;
	      ncoll++;
	    }
	    if (ncoll % 2 == 0) {
	      continue;
	    }
	    pixcoords[0] = px;
	    pixcoords[1] = py;
	  }
	  
	  //printf("pixcoords are %g %g %g\n",pixcoords[0],pixcoords[1],pixcoords[2]);

	  find_xi(ielem, pixcoords, xi_data, ei[pg->imtrx]->ielem_type, nodecoor, si, N_ext);
	  
	  pix_interp = 0;
	  for (ilnode = 0; ilnode < ei[pg->imtrx]->num_local_nodes; ilnode++) {
	    ignode = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ilnode];
	    phi_i = newshape(xi_data, ei[pg->imtrx]->ielem_type, PSI, ilnode, ei[pg->imtrx]->ielem_shape, efv->i[N_ext], ilnode);
	    pix_interp += phi_i*nodal_var_vals[ignode];
	  }
	  
	  tot_error += (pix_interp-pixdata[i][j][k])*(pix_interp-pixdata[i][j][k]);
	  numpoints++;
	}	
      }    
    }
  }
 
  tot_error /= numpoints;
  tot_error = sqrt(tot_error);
  return (tot_error);
}


static double bi_interp(double xg, double yg,
			double resolution[],
			double ***pixdata,
			int pixsize[])
{
/**************************************************************
Routine for bilinear interpolation (wikipedia has a nice entry explaining this)
Note that this has nothing to do with FEM, it's just a way to get the pixel value 
at an arbitrary position, by taking into account not just the pixel value at that position,
but also its nearest neighbors.

*/
  double dx, dy, dx2, dy2, dxdyinv;
  int ix, iy;
  double x1, x2, y1, y2;
  double xmax, ymax;
  double f11, f12, f21, f22;
  double val;

  dx = resolution[0];
  dy = resolution[1];
  dx2 = 0.5*dx;
  dy2 = 0.5*dy;
  dxdyinv = 1/(dx*dy);

  xmax = pixsize[0]*dx - dx2; 
  ymax = pixsize[1]*dy - dy2;

  ix = (int) ((xg - dx2)/ dx);
  iy = (int) ((yg - dy2)/ dy);
  
  if (xg < dx2) ix = -1;
  if (yg < dy2) iy = -1;

  x1 = MAX(dx2 + ix*dx,dx2);
  x2 = MIN(dx2 + (ix+1)*dx,xmax);
  y1 = MAX(dy2 + iy*dy,dy2);
  y2 = MIN(dy2 + (iy+1)*dy,ymax);
  f11 = pixdata[(int)(x1/dx)][(int)(y1/dy)][0];
  f12 = pixdata[(int)(x1/dx)][(int)(y2/dy)][0];
  f21 = pixdata[(int)(x2/dx)][(int)(y1/dy)][0];
  f22 = pixdata[(int)(x2/dx)][(int)(y2/dy)][0];
  val = dxdyinv*(f11*(x2-xg)*(y2-yg) + 
		 f21*(xg-x1)*(y2-yg) +
		 f12*(x2-xg)*(yg-y1) +
		 f22*(xg-x1)*(yg-y1));
    
  return val;
}
  

static double tri_interp(double xg, double yg, double zg,
			double resolution[],
			double ***pixdata,
			int pixsize[])
{
/**************************************************************
Routine for trilinear interpolation will go here. Same idea as above, but for 3D



  double dx, dy, dz;

  dx = resolution[0];
  dy = resolution[1];
  dz = resolution[2];
*/
  return 0;
}  

