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

/* ac_particles -- calculations involving discrete particles. 
 */
/*
 * $Id: ac_particles.c,v 5.6 2009-04-24 23:42:32 hkmoffa Exp $
 */

/* Needed to declare POSIX function drand48 */
#define _XOPEN_SOURCE

/* Standard include files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* GOMA include files */
#define GOMA_AC_PARTICLES_C
#include "goma.h"


/* Global variables extern declared in ac_particles.h. */
int Particle_Dynamics;		/* global toggle indicating particles are present. */
enum Particle_Model_t Particle_Model; /* What flavor of particle<->continuum stuff... */
dbl Particle_Model_Data[MAX_PARTICLE_MODEL_DATA_VALUES]; /* Real values for this model. */
int Particle_Number;		/* number of discrete particles. */
particle_filename_s Particle_Restart_Filename; /* restart filename */
int Particle_Output_Stride;	/* How often to output particle information. */
dbl Particle_Output_Time_Step;	/* Output every these units. */
int Particle_Max_Time_Steps;	/* max number of particle time steps if steady solution. */
enum Particle_Output_Format_t Particle_Output_Format; /* What kind of output file? */
int Particle_Full_Output_Stride;         /* > 0 => full output every that many steps. */
particle_filename_s Particle_Full_Output_Filename; /* where to put them. */
int Particle_Number_Sample_Types; /* How many datasets to output? */
int *Particle_Number_Samples_Existing; /* How many are tagged for each sample type?*/
int *Particle_Number_Samples;	/* How many particles to output for dataset #n? */
int *Particle_Number_Output_Variables; /* How many output vars for each sample. */
particle_variable_s **Particle_Output_Variables; /* List of variable indices to output for dataset #n */
particle_filename_s *Particle_Filename_Template; /* Template of where to put the data... */

dbl Particle_Density;		/* Density of particle in problem units */
dbl Particle_Radius;		/* Radius of particle in problem units. */
dbl Particle_Ratio;		/* Real/computational particle ratio. */
int Particle_Show_Debug_Info;	/* Show particle debug info. */
enum Particle_Domain_t Particle_Creation_Domain;
enum Particle_Domain_t Particle_Move_Domain;
particle_filename_s Particle_Creation_Domain_Filename;
particle_s Particle_Creation_Domain_Name;
particle_filename_s Particle_Move_Domain_Filename;
particle_s Particle_Move_Domain_Name;
dbl Particle_Creation_Domain_Reals[MAX_DOMAIN_REAL_VALUES];
dbl Particle_Move_Domain_Reals[MAX_DOMAIN_REAL_VALUES];
dbl xi_boundary_tolerances[3] = { XI_BOUNDARY_TOLERANCE0, XI_BOUNDARY_TOLERANCE1, XI_BOUNDARY_TOLERANCE2 };

int Particle_Number_PBCs;		/* Number of particle-related sideset BC's. */
PBC_t *PBCs;			/* Particle boundary condition structures. */

/* Global variables that reside entirely within this file. */
static particle_t *particles_to_do, *particles_to_send;
static particle_t **element_particle_list_head;
static particle_t **backup_element_particle_list_head;
static int num_particles;
static int last_element_loaded;
static int velo_interp;
static int last_elements_nodes_loaded;
static int last_elements_fv_loaded;
static dbl last_fvs_xi_loaded[DIM];
static my_fv_old_t *my_fv_old;

static int max_newton_iterations, max_particle_iterations;
static dbl total_accum_ust, output_accum_ust, particle_accum_ust;
#ifdef PARALLEL
static int total_num_particle_ghosts;
static dbl communication_accum_ust;
#endif

static int nodes_per_element;
static int sides_per_element;
static dbl **node_coord;

static const Exo_DB * static_exo;
static dbl * static_x;
static dbl * static_x_old;
static dbl * static_xdot;
static dbl * static_xdot_old;
static dbl * static_resid_vector;

element_particle_info_t *element_particle_info;

static dbl my_volume;
static dbl *el_volume;
static FILE **pa_fp;
static FILE *pa_full_fp;

static int pdim;		/* particle dimension (coordinate, velocities, etc.) */
static int mdim;		/* mesh dimension */

static int Num_Impermeable_PBCs;	/* Number of IMPERMEABLE particle boundary conditions. */

/* Function Prototypes */
static void find_exit_wound(const int, particle_t *, const dbl [3], const dbl [3], const dbl [3], const int, const int);

static dbl fill_element_volumes(const int, const int);

static int position_particle_uniformly(const int, particle_t *);

static char * construct_filename(const char *);

static void get_boundary_xi_newton(const dbl * const, const int, dbl *, const int, const dbl);

static particle_t * obtain_particle_space(const int);

static particle_t * create_a_particle(particle_t *, const int);

static int rejection_sample_a_particle(void);

static void initialize_and_create_new_particle(particle_t *, const int);

static void create_element_particle_info_maps(void);

static void handle_surface_interaction(particle_t *, dbl *, dbl *, int );

static void generate_source_particles(const dbl, const dbl, const dbl);

static void zero_a_particle(particle_t *);

static void initialize_surface_interactions(void);

static void load_element_information(const int);

static void load_element_node_coordinates(const int);

static void fill_hex_side_indices(const int, int[4]);

static void select_random_side_location(const int , const int, dbl *, dbl *);

static dbl get_side_area(const int, const int);

static void load_field_variables_at_xi(const int, dbl * const);

static int get_element_xi_newton(const int, const dbl * const, dbl *);

static dbl get_element_minimum_side_length(const int);

static void output_TECPLOT_zone_info(const dbl, const int, const int);

static void output_a_particle(particle_t * const, const dbl, const dbl, const int, const int, const int);

static void store_particle_data(particle_t * const);

static dbl compute_particle_dt(particle_t * const, const dbl);

static int move_particle(particle_t * const, dbl *, const dbl, const dbl);

static dbl my_distance(const dbl *, const dbl *, const int);

static void enforce_impermeable_PBCs(particle_t * const, dbl *, const dbl);

static dbl drand_standard_normal(void);

static dbl drand_truncated_normal(const dbl, const dbl);

static void dump_fcn(const char * const, const int, const int, particle_t *, ...);

static void fill_my_fv_old(void);

static void advance_a_particle(particle_t *, const dbl, const dbl, const int);

static void remove_from_element_particle_list(particle_t *, const int);

static void add_to_send_list(particle_t *);

static void add_to_do_list(particle_t *);

static void couple_to_continuum(void);

static void load_restart_file(void);

static void test_map_integrity(void);



/* Initialize the particles' positions and angular distributions.
 * Currently, only a uniform distribution is available.  This was a
 * pain in the ass to do...
 *
 * Unfortunately, this code assumes e_begin = 0... Oh well.
 */
void
initialize_particles(const Exo_DB * exo,
		     dbl * const x,
		     dbl * const x_old,
		     dbl * const xdot,
		     dbl * const xdot_old,
		     dbl * const resid_vector)
{
  dbl r, total_volume = 0.0;
  int  rejection, rejection_count, creation_count;
  int i, j, el_index_begin, el_index_end, num_els, elem_id;
  char time_of_day[80];
  particle_t p;

#ifdef PARALLEL
  int proc_index, left_over;
  int *number_to_create = 0, *local_number_created = 0, *number_created = 0;
  dbl *proc_volume_local = 0, *proc_volume = 0;
#endif

  DPRINTF(stdout, "\nParticle initialization:\n"); fflush(stdout);

  if((Particle_Model == SWIMMER_EXPLICIT || Particle_Model == SWIMMER_IMPLICIT) &&
     pd_glob[0]->CoordinateSystem != CARTESIAN)
    EH(-1, "You probably need to fix something first...");

  /* Allocate one of these.  I like the -> references instead of the
   * . references... */
  my_fv_old = (my_fv_old_t *)malloc(sizeof(my_fv_old_t));

  /* Set these at the two entries into particle stuff so I don't have
   * to pass them around EVERYWHERE. */
  static_exo = exo;
  static_x = x;
  static_x_old = x_old;
  static_xdot = xdot;
  static_xdot_old = xdot_old;
  static_resid_vector = resid_vector;

  /* Get the appropriate dimension for the coordinate, velocity loops... */
  pdim = pd_glob[0]->Num_Dim;
  if(pd_glob[0]->CoordinateSystem == SWIRLING ||
     pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN ||
     pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D)
    pdim = pdim + 1;
  
  /* Get mesh dimension */
  mdim = static_exo->num_dim;

  if(pd_glob[0]->e[pg->imtrx][R_MESH1])
    EH(-1, "Cannot couple particles and deformable mesh, yet.");
  if(mdim == 2)
    {
      if(static_exo->eb_elem_itype[0] != BIQUAD_QUAD &&
	 static_exo->eb_elem_itype[0] != BILINEAR_QUAD)
	EH(-1, "Can only have 9-node and 4-node quadrilateral elements in 2D.");
      sides_per_element = 4;
    }
  else /* static_exo->num_dim == mdim == 3 */
    {
      if(static_exo->eb_elem_itype[0] != TRIQUAD_HEX &&
	 static_exo->eb_elem_itype[0] != TRILINEAR_HEX)
	EH(-1, "Can only handly 27-node and 8-node hex elements in 3D.");
      sides_per_element = 6;
    }

  nodes_per_element = elem_info(NNODES, static_exo->eb_elem_itype[0]);
  node_coord = (dbl **)malloc(nodes_per_element * sizeof(dbl *));
  for(i = 0; i < nodes_per_element; i++)
    node_coord[i] = (dbl *)malloc(DIM * sizeof(dbl));

  element_particle_list_head = (particle_t **)malloc(static_exo->num_elems * sizeof(particle_t *));
  backup_element_particle_list_head = (particle_t **)malloc(static_exo->num_elems * sizeof(particle_t *));
  particles_to_send = NULL;

  DPRINTF(stdout, "  elements: "); fflush(stdout);
  /* We will always allocate this element-sized space even though we
   * might not always use it.  There are too many use cases (some
   * delayed) that complicate this... */
  element_particle_info = (element_particle_info_t *)malloc(static_exo->num_elems * sizeof(element_particle_info_t));
  for(i = 0; i < static_exo->num_elems; i++)
    {
      element_particle_info[i].PBC_side_id = (int *)malloc(sides_per_element * sizeof(int));
#ifdef PARALLEL
      element_particle_info[i].owner_local_element_id = (int *)calloc((unsigned)sides_per_element, sizeof(int));
      element_particle_info[i].num_ghost_target_elems = 0; /* number of times each particle needs to be ghosted. */
      element_particle_info[i].ghost_proc = NULL; /* what proc do they ghost to? */
      element_particle_info[i].ghost_local_elem_id = NULL; /* what's the local element id there? */
#endif
      for(j = 0; j < sides_per_element; j++)
	{
	  element_particle_info[i].PBC_side_id[j] = -1;
#ifdef PARALLEL
	  element_particle_info[i].owner_local_element_id[j] = -1;
#endif
	}
      if(Particle_Model == SWIMMER_IMPLICIT || Particle_Model == SWIMMER_EXPLICIT)
	element_particle_info[i].source_term = (dbl *)calloc(MDE, sizeof(dbl));
      else
	element_particle_info[i].source_term = NULL;
      element_particle_info[i].list_PBC_IMPERMEABLE = NULL;
      element_particle_list_head[i] = NULL;
      backup_element_particle_list_head[i] = NULL;
    }
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

  DPRINTF(stdout, "  output: "); fflush(stdout);
  /* Open the output file for full output/restart if we've selected one. */
  if(Particle_Full_Output_Stride)
    {
      if(!(pa_full_fp = fopen(construct_filename(Particle_Full_Output_Filename), "w")))
	EH(-1, "Could not open Particle_Full_Output_Filename for zeroing.");
      DPRINTF(stdout, "%s, ", Particle_Full_Output_Filename); fflush(stdout);
    }
  else
    pa_full_fp = NULL;

  if(Particle_Number_Sample_Types)
    {
      pa_fp = (FILE **)array_alloc(1, Particle_Number_Sample_Types, sizeof(FILE *));
      for(i = 0; i < Particle_Number_Sample_Types; i++)
	{
	  if(!(pa_fp[i] = fopen(construct_filename(Particle_Filename_Template[i]), "w")))
	    EH(-1, "Could not open a particle file for zeroing.");
	  DPRINTF(stdout, "%s, ", Particle_Filename_Template[i]); fflush(stdout);
	  Particle_Number_Samples_Existing[i] = 0;
	  if(Particle_Output_Format == TECPLOT)
	    {
	      fprintf(pa_fp[i], "TITLE = \"Goma Particles\"\n");
#ifdef PARALLEL
	      if(pdim == 2)
		fprintf(pa_fp[i], "VARIABLES = \"X\", \"Y\", \"PROCID\", \"TIME\"");
	      else
		fprintf(pa_fp[i], "VARIABLES = \"X\", \"Y\", \"Z\", \"PROCID\", \"TIME\"");
#else
	      if(pdim == 2)
		fprintf(pa_fp[i], "VARIABLES = \"X\", \"Y\", \"TIME\"");
	      else
		fprintf(pa_fp[i], "VARIABLES = \"X\", \"Y\", \"Z\", \"TIME\"");
#endif
	      for(j = 0; j < Particle_Number_Output_Variables[i]; j++)
		fprintf(pa_fp[i], ", \"%s\"", Particle_Output_Variables[i][j]);
	      fprintf(pa_fp[i], "\n");
	    }
	}
    }
  else
    pa_fp = NULL;
  if(Particle_Output_Format == TECPLOT &&
     Particle_Number > 0)
    output_TECPLOT_zone_info(0.0, 0, 1);
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

  last_element_loaded = -1;
  last_elements_nodes_loaded = -1;
  last_elements_fv_loaded = -1;
  for(i = 0; i < DIM; i++)
    last_fvs_xi_loaded[i] = -1.0e+10;

  num_particles = 0;

  DPRINTF(stdout, "  volumes: "); fflush(stdout);
  /* Position initial particles. */
  if(Particle_Number > 0)
    {
      el_index_begin = static_exo->eb_ptr[0];
      el_index_end = static_exo->eb_ptr[static_exo->num_elem_blocks];
      num_els = el_index_end - el_index_begin;
      el_volume = (dbl *)malloc(num_els * sizeof(dbl));
  
      /* First get total element "volume" */
      total_volume = my_volume =  fill_element_volumes(el_index_begin, el_index_end);
#ifdef PARALLEL
      proc_volume = (dbl *)calloc((unsigned)Num_Proc, sizeof(dbl));
      proc_volume_local = (dbl *)calloc((unsigned)Num_Proc, sizeof(dbl));
      proc_volume_local[ProcID] = my_volume;
      MPI_Reduce(proc_volume_local, proc_volume, Num_Proc, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
      /* Get total volume. */
      total_volume = 0.0;
      for(i = 0; i < Num_Proc; i++)
	total_volume += proc_volume[i];

      /*
	DPRINTF(stderr, "Total volume = %g\n", total_volume); fflush(stderr);
      */
#endif
    }
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

  DPRINTF(stdout, "  domains: "); fflush(stdout);
  if(Particle_Number > 0)
    {
      /* This is really the only check we can make unless I want to
       * implement some pretty serious probabilistic sampling and
       * in/out comparisons for ACIS objects. */
      if(Particle_Creation_Domain == BRICK && Particle_Move_Domain == BRICK)
	{
	  for(i = 0; i < 3; i++)
	    if(Particle_Move_Domain_Reals[2*i] > Particle_Creation_Domain_Reals[2*i] ||
	       Particle_Move_Domain_Reals[2*i+1] < Particle_Creation_Domain_Reals[2*i+1])
	      EH(-1, "Mismatch in BRICK bounds for Creation/Move.");
	}
    }
  
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

  DPRINTF(stdout, "  particles: "); fflush(stdout);

  /* If Particle_Number == -1, then we're performing the "particle
   * mass test".  One particle is placed on each element, and the
   * total particle mass is interpolated by shape functions.  The
   * global value should sum to be, of course, (# particles) x (# real
   * particles per computational particle) x (mass per real particle).
   * Note that in this situation # particles = # elements. */
  if(Particle_Number == -1)
    {
      el_index_begin = static_exo->eb_ptr[0];
      el_index_end = static_exo->eb_ptr[static_exo->num_elem_blocks];
      num_els = el_index_end - el_index_begin;
      for(elem_id = el_index_begin; elem_id < el_index_end; elem_id++)
	{
#ifdef PARALLEL
	  if(DPI_ptr->elem_owner[elem_id] != ProcID)
	    continue;
#endif

	  /* Place the particles at the arithmetic mean of the node
	   * coordinates.  This will place the particle in the middle
	   * of the element as long as the element is non-concave. */
	  load_element_node_coordinates(elem_id);
	  zero_a_particle(&p);
	  for(i = 0; i < mdim; i++)
	    {
	      for(j = 0; j < nodes_per_element; j++)
		p.x[i] += node_coord[j][i];
	      p.x[i] /= (dbl)nodes_per_element;
	    }
	  
	  initialize_and_create_new_particle(&p, elem_id);
	}
    }

  if(Particle_Number > 0)
    {
      rejection_count = 0;
      creation_count = 0;
#ifdef PARALLEL
      DPRINTF(stdout, "\n"); fflush(stdout);
      /* Processor zero controls the global random sampling procedure.
       * This is pretty trivial in single-processor. */
      number_to_create = (int *)calloc((unsigned)Num_Proc, sizeof(int));
      number_created = (int *)calloc((unsigned)Num_Proc, sizeof(int));
      local_number_created = (int *)calloc((unsigned)Num_Proc, sizeof(int));
      if(ProcID == 0)
	{
	  while(creation_count < Particle_Number)
	    {
	      DPRINTF(stdout, "    %d=", Particle_Number - creation_count); fflush(stdout);
	      /* First figure out how much each processor should generate, rounded down. */
	      left_over = Particle_Number - creation_count;
	      for(proc_index = 0; proc_index < Num_Proc; proc_index++)
		{
		  number_to_create[proc_index] = (int)floor(proc_volume[proc_index]/total_volume * (Particle_Number - creation_count));
		  left_over -= number_to_create[proc_index];
		}

	      /* Now we need to assign the leftovers. */
	      for(i = 0; i < left_over; i++)
		{
		  /* This random number determines which processor we're going
		   * to try to create a particle on. */
		  r = drand48() * total_volume;
		  proc_index = 0;
		  while(proc_index < Num_Proc &&
			r > proc_volume[proc_index])
		    {
		      r -= proc_volume[proc_index];
		      proc_index++;
		    }
		  if(proc_index == Num_Proc)
		    EH(-1, "proc_index == Num_Proc");
		  number_to_create[proc_index]++;
		}

	      for(i = 0; i < Num_Proc - 1; i++)
		DPRINTF(stdout, "%d+", number_to_create[i]);
	      DPRINTF(stdout, "%d", number_to_create[Num_Proc - 1]); fflush(stdout);
		
	      /* Now broadcast the vector containing how many particles each processor should create. */
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Bcast(number_to_create, Num_Proc, MPI_INT, 0, MPI_COMM_WORLD);

	      /* Create our particles while everyone else is creating theirs. */
	      rejection = 0;
	      for(i = 0; i < number_to_create[0]; i++)
		rejection += rejection_sample_a_particle();

	      /* Now receive how many particles were created by everyone. */
	      memset(local_number_created, 0, Num_Proc * sizeof(int));
	      local_number_created[0] = number_to_create[0] - rejection;
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Reduce(local_number_created, number_created, Num_Proc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	      for(i = 0; i < Num_Proc; i++)
		{
		  creation_count += number_created[i];
		  rejection_count += number_to_create[i] - number_created[i];
		}
	      DPRINTF(stdout, " -> ");
	      for(i = 0; i < Num_Proc - 1; i++)
		DPRINTF(stdout, "%d+", number_created[i]);
	      DPRINTF(stdout, "%d", number_created[Num_Proc - 1]); fflush(stderr);
	      DPRINTF(stdout, ", total C/R=%d/%d\n", creation_count, rejection_count);
	    }
	  number_to_create[0] = -1;
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Bcast(number_to_create, Num_Proc, MPI_INT, 0, MPI_COMM_WORLD);
	}
      else
	{
	  while(number_to_create[0] != -1)
	    {
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Bcast(number_to_create, Num_Proc, MPI_INT, 0, MPI_COMM_WORLD);
	      if(number_to_create[0] == -1)
		continue;
	      rejection = 0;
	      for(i = 0; i < number_to_create[ProcID]; i++)
		rejection += rejection_sample_a_particle();
	      memset(local_number_created, 0, Num_Proc * sizeof(int));
	      local_number_created[ProcID] = number_to_create[ProcID] - rejection;
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Reduce(local_number_created, number_created, Num_Proc, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	    }
	}
      DPRINTF(stdout, "    total R = %d, ", rejection_count); fflush(stdout);
#else
      /* Single processor version. */      
      while(creation_count < Particle_Number)
	{
	  rejection = rejection_sample_a_particle();
	  rejection_count += rejection;
	  creation_count += (1 - rejection);
	}
      DPRINTF(stdout, "total R=%d, ", creation_count); fflush(stdout);
#endif
    }
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

  DPRINTF(stdout, "  restart: "); fflush(stdout);
  if(strncmp(Particle_Restart_Filename, "<not active>", 12))
    {
      load_restart_file();
      DPRINTF(stdout, "%s, ", Particle_Restart_Filename);
    }
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);
  
  DPRINTF(stdout, "  surfaces: "); fflush(stdout);
  initialize_surface_interactions();
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);
  
  DPRINTF(stdout, "  maps: "); fflush(stdout);
  create_element_particle_info_maps();
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);
  
  DPRINTF(stdout, "  coupling: "); fflush(stdout);
  couple_to_continuum();
#ifdef PARALLEL
  DPRINTF(stdout, "%d ghosts, ", total_num_particle_ghosts);
#endif
  get_time(time_of_day);
  DPRINTF(stdout, "done at %s\n", time_of_day); fflush(stdout);

#ifdef PARALLEL
  if(Particle_Number > 0)
    {
      free(proc_volume);
      free(proc_volume_local);
      free(number_to_create);
      free(number_created);
      free(local_number_created);
    }
#endif
}


/* This routine should be called to find out what side a particle
 * intersected as it left the problem domain.  It handles cases where
 * the particle intersected exactly one extended side on its way out,
 * as well as the case of the particle intersecting two extended sides
 * on its way out.  In either case, it computes the point of
 * intersection on the boundary by using Newton's method.
 *
 * All calls to get_boundary_xi_newton also modify fv->x, etc.
 *
 * Upon entering, p is a pointer to the particle in question, and
 * old_el_id is the element it was in before it exited the domain.
 *
 * This is where the final particle location is computed, especially
 * with respect to particle boundary conditions. */
static void
find_exit_wound(const int elem_id,
		particle_t * p,
		const dbl x_start[DIM],
		const dbl x_end[DIM],
		const dbl xi[DIM],
		const int stack_count,
		const int tolerance_level)
{
  dbl xi_tmp[DIM], x_intersect[DIM], first_xi_tmp[DIM], second_xi_tmp[DIM] = {0.0, 0.0, 0.0};
  dbl coeff[6];
  /* dbl len; */
  int exit_id1, exit_id2, exit_id3, PBC_id;
  int bdry_crossed[3], num_crossed;
  int i;
  /* initialize first_xi_tmp */
  for (i = 0; i < DIM; i++) {
    first_xi_tmp[i] = 0;
  }

  if(stack_count == 100)
    dump1(EXIT, p, "Hit 100 recursive calls to find_exit_wound(), something must be wrong...");

  exit_id1 = exit_id2 = exit_id3 = -1;

  /* Not necessary except for debugging. */
  load_element_node_coordinates(elem_id);

  if(mdim == 2)
    {
      /* For 2D, these are always 0.0. */
      xi_tmp[2] = x_intersect[2] = 0.0;
      memcpy(xi_tmp, xi, 2 * sizeof(dbl));

      /* Write the physical coordinate line segment as a line, ax + by + c = 0 */
      coeff[0] = x_start[1] - x_end[1];
      coeff[1] = x_end[0] - x_start[0];
      coeff[2] = -(coeff[0] * x_start[0] + coeff[1] * x_start[1]);
      
      bdry_crossed[0] = (fabs(xi[0]) > 1.0) ? 1 : 0;
      bdry_crossed[1] = (fabs(xi[1]) > 1.0) ? 1 : 0;
      bdry_crossed[2] = 0;
      num_crossed = bdry_crossed[0] + bdry_crossed[1];

      if(num_crossed == 2)
	{
	  /* Find the two intersecting sides (0-based). */
	  exit_id1 = (xi[0] < -1.0) ? 3 : 1;
	  exit_id2 = (xi[1] < -1.0) ? 0 : 2;
      
	  /* First, check the intersection corresponding to exit_id1. */
	  xi_tmp[0] = (exit_id1 == 3 ? -1.0 : 1.0);
	  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, 1.0e-12);

	  if(fabs(xi_tmp[1]) > 1.0) /* Not this intersection!  Do exit_id2... */
	    {
	      xi_tmp[0] = xi[0]; /* reset initial guess for Newton's method. */
	      xi_tmp[1] = (exit_id2 == 0 ? -1.0 : 1.0);
	      get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, 1.0e-12);

	      if(fabs(xi_tmp[0]) > 1.0)
		dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
		       x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
		       bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
		       xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);

	      exit_id1 = exit_id2; /* so exit_side has the right value. */
	    }
	}
      else
	{
	  /* At this point only one of xi_tmp[*] is < -1.0 or > 1.0. */
	  if(bdry_crossed[0])
	    {
	      if(xi[0] > 1.0) /* to side 1 (0-based) */
		{
		  xi_tmp[0] = 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, 1.0e-12);
		  exit_id1 = 1;
		}
	      else /* xi[0] < -1.0, to side 3 (0-based) */
		{
		  xi_tmp[0] = -1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, 1.0e-12);
		  exit_id1 = 3;
		}
	    }
	  else /* fabs(xi[1]) > 1.0 */
	    {
	      if(xi[1] > 1.0) /* to side 2 (0-based) */
		{
		  xi_tmp[1] = 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, 1.0e-12);
		  exit_id1 = 2;
		}
	      else /* xi[1] < -1.0, to side 0 (0-based) */
		{
		  xi_tmp[1] = -1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, 1.0e-12);
		  exit_id1 = 0;
		}
	    }
	}
    }
  else
    {  /* 3D crossings */
      /* Write the physical coordinate line segment as a line, x_start + t * (x_end - x_start) = 0 */
      for(i = 0; i < 3; i++)
	{
	  coeff[i] = x_start[i];
	  coeff[i+3] = x_end[i] - x_start[i];
	}
      /* len = nnorm(3, &(coeff[3])); */

      /* We know len > 0.0 b/c we crossed a boundary to get into this
       * routine; hence, we must have moved. */
      /*
      for(i = 3; i < 6; i++)
	coeff[i] /= len;
      */
      
      /* Initialize ... */
      memcpy(xi_tmp, xi, 3 * sizeof(dbl));

      /* Calculate which xi[0,1,2] faces we crossed */
      num_crossed = 0;
      for(i = 0; i < 3; i++)
	{
	  bdry_crossed[i] = (fabs(xi[i]) > 1.0) ? 1 : 0;
	  num_crossed += bdry_crossed[i];
	}
      
      /*
      fprintf(stderr, "Crossed %d boundaries: %d %d %d\n", num_crossed,
	      bdry_crossed[0], bdry_crossed[1], bdry_crossed[2]);
      fprintf(stderr, "  xi = (%g,%g,%g)\n", xi[0], xi[1], xi[2]);
      dump_fcn(__FILE__, __LINE__, NO_EXIT, p, "");
      */
      
      if(num_crossed == 3)
	{
	  /* We are in one of the eight "triple" crossings. */

	  /* Find the three intersecting faces (0-based). */
	  exit_id1 = (xi[0] < -1.0) ? 3 : 1;
	  exit_id2 = (xi[1] < -1.0) ? 0 : 2; /* possibly needs to be switched */
	  exit_id3 = (xi[2] < -1.0) ? 4 : 5;
      
	  /* First, check the intersection corresponding to exit_id1. */
	  xi_tmp[0] = (exit_id1 == 3 ? -1.0 : 1.0);
	  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, xi_boundary_tolerances[tolerance_level]);

	  if(fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0) /* Not this intersection!  Do exit_id2... */
	    {
	      xi_tmp[0] = xi[0]; /* reset initial guess for Newton's method. */
	      xi_tmp[2] = xi[2];
	      xi_tmp[1] = (exit_id2 == 0 ? -1.0 : 1.0);
	      get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, xi_boundary_tolerances[tolerance_level]);

	      if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[2]) > 1.0) /* Ok,third one's a charm, right? */
		{
		  xi_tmp[0] = xi[0]; /* reset initial guess for Newton's method. */
		  xi_tmp[1] = xi[1];
		  xi_tmp[2] = (exit_id3 == 4 ? -1.0 : 1.0);
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 2, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  fprintf(stderr, "second xi_tmp = (%g,%g,%g)\n", second_xi_tmp[0], second_xi_tmp[1], second_xi_tmp[2]);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Trying again...\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  fprintf(stderr, "second xi_tmp = (%g,%g,%g)\n", second_xi_tmp[0], second_xi_tmp[1], second_xi_tmp[2]);
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = exit_id3; /* exit_id1 should have the intersecting face. */
		}
	      else
		exit_id1 = exit_id2; /* exit_id1 should have the final intersecting face. */
	    }
	}
      else if(num_crossed == 2)
	{
	  if(!bdry_crossed[0]) /* crossed xi[1] and xi[2] */
	    {
	      exit_id1 = (xi[1] < -1.0) ? 0 : 2;
	      exit_id2 = (xi[2] < -1.0) ? 4 : 5;

	      /* Check exit_id1 first. */
	      xi_tmp[1] = (exit_id1 == 0) ? -1.0 : 1.0;
	      get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, xi_boundary_tolerances[tolerance_level]);
	      /*
	      fprintf(stderr, "xi_tmp = (%g,%g,%g)\n", xi_tmp[0], xi_tmp[1], xi_tmp[2]);
	      */
	      if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[2]) > 1.0) /* Not this intersction!  Do exit_id2... */
		{
		  memcpy(first_xi_tmp, xi_tmp, DIM*sizeof(dbl));
		  xi_tmp[1] = xi[1]; /* reset initial guess for Newton's method. */
		  xi_tmp[2] = (exit_id2 == 4) ? -1.0 : 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 2, xi_boundary_tolerances[tolerance_level]);

		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Trying again...\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  else
		    exit_id1 = exit_id2; /* exit_id1 should have the intersecting face. */
		}
	    }
	  else if(!bdry_crossed[1]) /* crossed xi[0] and xi[2] */
	    {
	      exit_id1 = (xi[0] < -1.0) ? 3 : 1;
	      exit_id2 = (xi[2] < -1.0) ? 4 : 5;

	      /* Check exit_id1 first. */
	      xi_tmp[0] = (exit_id1 == 3) ? -1.0 : 1.0;
	      get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, xi_boundary_tolerances[tolerance_level]);
	      
	      if(fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0) /* Not this intersction!  Do exit_id2... */
		{
		  memcpy(first_xi_tmp, xi_tmp, DIM*sizeof(dbl));
		  xi_tmp[0] = xi[0]; /* reset initial guess for Newton's method. */
		  xi_tmp[2] = (exit_id2 == 4) ? -1.0 : 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 2, xi_boundary_tolerances[tolerance_level]);

		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Trying again...\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  else
		    exit_id1 = exit_id2; /* exit_id1 should have the intersecting face. */
		}
	    }
	  else /* crossed xi[0] and xi[1] */
	    {
	      exit_id1 = (xi[0] < -1.0) ? 3 : 1;
	      exit_id2 = (xi[1] < -1.0) ? 0 : 2;

	      /* Check exit_id1 first. */
	      xi_tmp[0] = (exit_id1 == 3) ? -1.0 : 1.0;
	      get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, xi_boundary_tolerances[tolerance_level]);
	      
	      if(fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0) /* Not this intersction!  Do exit_id2... */
		{
		  memcpy(first_xi_tmp, xi_tmp, DIM*sizeof(dbl));
		  xi_tmp[0] = xi[0]; /* reset initial guess for Newton's method. */
		  xi_tmp[1] = (exit_id2 == 0) ? -1.0 : 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, xi_boundary_tolerances[tolerance_level]);

		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Trying again...\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  fprintf(stderr, "first xi_tmp = (%g,%g,%g)\n", first_xi_tmp[0], first_xi_tmp[1], first_xi_tmp[2]);
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  else
		    exit_id1 = exit_id2; /* exit_id1 should have the intersecting face. */
		}
	    }
	}
      else /* one coordinate crossing. */
	{
	  if(bdry_crossed[0]) /* crossed a xi[0] face */
	    {
	      if(xi[0] <= -1.0) /* xi[0] = -1 */
		{
		  xi_tmp[0] = -1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Trying again...\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 3;
		}
	      else /* xi[0] = +1 */
		{
		  xi_tmp[0] = 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 0, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 1;
		}
	    }
	  else if(bdry_crossed[1]) /* crossed a xi[1] face */
	    {
	      if(xi[1] <= -1.0) /* xi[1] = -1 */
		{
		  xi_tmp[1] = -1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 0;
		}
	      else /* xi[1] = +1 */
		{
		  xi_tmp[1] = 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 1, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 2;
		}
	    }
	  else /* crossed a xi[2] face */
	    {
	      if(xi[2] <= -1.0) /* xi[2] = -1 */
		{
		  xi_tmp[2] = -1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 2, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 4;
		}
	      else /* xi[2] = +1 */
		{
		  xi_tmp[2] = 1.0;
		  get_boundary_xi_newton(coeff, elem_id, xi_tmp, 2, xi_boundary_tolerances[tolerance_level]);
		  if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0)
		    {
		      if(tolerance_level < 2)
			{
			  if(Num_Proc > 1)
			    fprintf(stderr, "WARNING: Proc%d: Going to tolerance level %d in find_exit_wound().\n",
				    ProcID, tolerance_level + 1);
			  else
			    fprintf(stderr, "WARNING: Going to tolerance level %d in find_exit_wound().\n",
				    tolerance_level + 1);
			  dump20(NO_EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			  find_exit_wound(elem_id, p, x_start, x_end, xi, stack_count+1, tolerance_level+1);
			  return;
			}
		      else
			{
			  dump20(EXIT, p, "  x_start = (%g,%g,%g), x_end = (%g,%g,%g)\n  bdry_crossed = (%d,%d,%d), exit_ids = (%d,%d,%d)\n  xi = (%g,%g,%g), xi_tmp = (%g,%g,%g)\n  stack_count = %d\nOh no!  I couldn't find any good intersections!  Bye-bye!\n",
				 x_start[0], x_start[1], x_start[2], x_end[0], x_end[1], x_end[2],
				 bdry_crossed[0], bdry_crossed[1], bdry_crossed[2], exit_id1, exit_id2, exit_id3,
				 xi[0], xi[1], xi[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], stack_count);
			}
		    }
		  exit_id1 = 5;
		}
	    }
	}
    }
  
  PBC_id = element_particle_info[elem_id].PBC_side_id[exit_id1];
  
  memcpy(x_intersect, fv->x, DIM * sizeof(dbl));

  /*
  fprintf(stderr, "  Intersection occured at x = (%g,%g,%g) and xi = (%g,%g,%g)\n",
	  fv->x[0], fv->x[1], fv->x[2], xi_tmp[0], xi_tmp[1], xi_tmp[2]);
  fprintf(stderr, "  PBC index for this intersection is %d (elem = %d, side = %d)\n",
	  PBC_id, elem_id, exit_id1);
  */
  
  if(PBC_id == -1)		/* must be another local element on the other side! */
    {
      /* Most of these PBC_id == -1 cases are merely particles moving
       * from one element to another.  If that's the case, update the
       * particle's element, and recalculate xi according to that
       * element.  If we're still out-of-element, then hit the next
       * part. */
      p->owning_elem_id = static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + exit_id1];

      if(p->owning_elem_id == -1)
	{
	  fprintf(stderr, "\nelem_id = %d, global_elem_id = %d\n", elem_id, DPI_ptr->elem_index_global[elem_id]);
	  fprintf(stderr, "elem_elem_list for this element is: %d %d %d %d %d %d\n",
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 0],
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 1],
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 2],
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 3],
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 4],
		  static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + 5]);
	  dump3(EXIT, p, "elem_id = %d, exit_id1 = %d\nEnded up in element -1!", elem_id, exit_id1);
	}

      if(get_element_xi_newton(p->owning_elem_id, x_end, xi_tmp) == -1)
	dump10(EXIT, p, "  exit_ids = (%d,%d,%d)\n  bdry_crossed = (%d,%d,%d)\n  xi_tmp = (%g,%g,%g)\nOh no!  I couldn't find element coordinates!  Bye-bye!\n",
	     exit_id1, exit_id2, exit_id3,
	     bdry_crossed[0], bdry_crossed[1], bdry_crossed[2],
	     xi_tmp[0], xi_tmp[1], xi_tmp[2]);

#ifdef SHOW_PARTICLE_MOVEMENT
      fprintf(stderr, "%d(%d)->", elem_id, exit_id1);
#endif

      /*
	fprintf(stderr, "  Side %d of element %d has no particle boundary condition (PBC)!\n",
	exit_id1, elem_id);
	fprintf(stderr, "  Checking (recursively) from element %d.\n", p->owning_elem_id);
      */

      /* Now this gets a little tricky, er, recursive.  Since we've
       * apparently jumped out of the domain by skipping over an
       * element (or possibly more than one), we pretend that we're
       * really in the adjacent element, with the same exit location,
       * and call this routine again.  It will make a recursive call
       * for each element we've "jumped" until we acutally look from
       * the boundary element, in which case this routine should
       * terminate just fine.  I could use this same logic to catch
       * "corner clippings" if I like.  Just something to keep in
       * mind... */
      if(fabs(xi_tmp[0]) > 1.0 || fabs(xi_tmp[1]) > 1.0 || fabs(xi_tmp[2]) > 1.0)
	{
	  last_element_loaded = -1;
	  last_elements_fv_loaded = -1;
	  find_exit_wound(p->owning_elem_id, p, x_start, x_end, xi_tmp, stack_count + 1, 0);
	  /*
	    find_exit_wound(p->owning_elem_id, p, x_intersect, x_end, xi_tmp);
	  */
	}
      else
	{
#ifdef SHOW_PARTICLE_MOVEMENT
	  fprintf(stderr, "%d\n", p->owning_elem_id);
#endif
          if (p->x != x_end) {
            memcpy(p->x, x_end, DIM * sizeof(dbl));
          }
	  memcpy(p->xi, xi_tmp, DIM * sizeof(dbl));
	}
    }
#ifdef PARALLEL    
  else if(PBC_id < -1)
      {
	/* Pack and ship this particle to another processor */
	memcpy(p->x_start, x_start, DIM * sizeof(dbl));
	memcpy(p->x_end, x_end, DIM * sizeof(dbl));
	p->owning_elem_id = element_particle_info[elem_id].owner_local_element_id[exit_id1];
	p->owning_proc_id = -(PBC_id+2);
	p->state = PROC_TRANSFER;
#ifdef SHOW_PARTICLE_MOVEMENT
	fprintf(stderr, "%d(%d) to P%d/%d\n", elem_id, exit_id1, -(PBC_id+2), p->owning_elem_id);
#endif
      }
#endif
  else
    {
      handle_surface_interaction(p, x_intersect, xi_tmp, PBC_id);
#ifdef SHOW_PARTICLE_MOVEMENT
      fprintf(stderr, "%d(%d) PBC/%d\n", elem_id, exit_id1, PBC_id);
#endif
    }
}


/* This function is equivalent in a general way to the assemble_*
 * functions for the regular FEM variables.  This routine moves the
 * particles by advection with a simple Weiner process for random
 * components.
 */
int
compute_particles(const Exo_DB *exo,	/* mesh */
		  dbl * const x,	/* solution vector */
		  dbl * const x_old,	/* old solution vector */
		  dbl * const xdot,	/* solution time derivative */
		  dbl * const xdot_old, /* old time derivative */
		  dbl * const resid_vector, /* residual vector */
		  const dbl global_end_time,	/* end of Goma's timestep */
                  const dbl goma_dt,	/* time step size */
		  const int n)	/* time step number */
{
  double start_ust, temp_ust, global_start_time;
  particle_t *p, *p_tmp;
  int i, el_index;
#ifdef PARALLEL
  int mpi_byte_count, mpi_retval;
  particle_t termination_p, send_p, recv_p;
  int j, done, local_max_particle_iterations, local_max_newton_iterations;
  int local_num_particles, local_particle_transfers, particle_transfers;
  dbl local_total_accum_ust, local_particle_accum_ust, local_output_accum_ust, local_communication_accum_ust;
  int local_num_to_send, num_to_send;
  MPI_Status mpi_status;
#endif
  
  total_accum_ust = 0.0;
  particle_accum_ust = 0.0;
  output_accum_ust = 0.0;
#ifdef PARALLEL
  local_particle_transfers = 0;
  communication_accum_ust = 0.0;
  num_to_send = 0;
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  start_ust = ust();
  
  /* Set these at the two entries into particle stuff so I don't have
   * to pass them around EVERYWHERE. */
  static_exo = exo;
  static_x = x;
  static_x_old = x_old;
  static_xdot = xdot;
  static_xdot_old = xdot_old;
  static_resid_vector = resid_vector;

  if(Particle_Output_Format == TECPLOT)
    output_TECPLOT_zone_info(global_end_time, n, 0);
  output_accum_ust = MAX(ust() - start_ust, 0.0);

  /* Initialize where we most recently looked up values. */
  last_element_loaded = -1;
  last_elements_nodes_loaded = -1;
  last_elements_fv_loaded = -1;
  for(i = 0; i < DIM; i++)
    last_fvs_xi_loaded[i] = 1.0e+10;

  /* Model-specific initializations. */
  if(Particle_Model == SWIMMER_EXPLICIT ||
     Particle_Model == SWIMMER_IMPLICIT)
    for(i = 0; i < static_exo->num_elems; i++)
      memset(element_particle_info[i].source_term, 0, MDE * sizeof(dbl));

  global_start_time = global_end_time - goma_dt;
  max_particle_iterations = 0;
  max_newton_iterations = 0;
#ifdef PARALLEL
  local_num_to_send = 0;
#endif
  
  /* Doing this at the begining of the time step.  Could be done at
   * the end, too. */
  generate_source_particles(global_end_time, goma_dt, global_start_time);

  /* First move all the particles that we already have. */
  for(el_index = 0; el_index < static_exo->num_elems; el_index++)
    {
#ifdef PARALLEL
      /* This shouldn't affect anything unless I was holding on to
       * ghosted particles somehow... */
      if(DPI_ptr->elem_owner[el_index] != ProcID)
	continue;
#endif
      
      p = element_particle_list_head[el_index];

      while(p)
	{
	  if(get_element_xi_newton(el_index, p->x, p->xi) == -1)
	    dump1(EXIT, p, "Sanity check.");
	  
	  if(Particle_Show_Debug_Info)
	    {
	      if(fabs(p->xi[0]) > 1.0 || fabs(p->xi[1]) > 1.0 || fabs(p->xi[2]) > 1.0)
		dump1(NO_EXIT, p, "PARTICLE NOT LOCATED ON PROPER ELEMENT");

	      if(p->owning_elem_id != el_index)
		dump1(NO_EXIT, p, "PARTICLE NOT LOCATED ON ITS OWNING ELEMENT");
	    }

	  p_tmp = p->next;	/* Save in case p gets deleted */
	  /*
	  p->time = global_start_time;
	  */
	  advance_a_particle(p, global_start_time, global_end_time, n);

	  if(p->state == PROC_TRANSFER) /* I need to go to another processor. */
	    {
	      remove_from_element_particle_list(p, el_index);
	      add_to_send_list(p);
	      num_particles--;
#ifdef PARALLEL
	      local_num_to_send++;
#endif
	      /*
	      fprintf(stderr, "Proc %d is adding a processor to its send list\n", ProcID);
	      fflush(stderr);
	      */
	    }
	  else if(p->state == DEAD) /* Delete me */
	    {
	      num_particles--;
	      remove_from_element_particle_list(p, el_index);
	      free(p);	/* Bye bye */
	      fprintf(stderr, "REMOVING A PARTICLE\n");
	    }
	  else if(p->state == ACTIVE && el_index != p->owning_elem_id) /* I moved elements on the same processor. */
	    {
	      /* This particle won't get redone because it's p->time is used. */
	      remove_from_element_particle_list(p, el_index);
	      create_a_particle(p, p->owning_elem_id);
	      num_particles--; /* This isn't really a new particle... */
	    }

	  if(p->state == ACTIVE)
	    {
	      /*
	      output_start_ust = ust();
	      if(TimeIntegration == STEADY)
		{
		  if(p->time == global_end_time)
		    output_a_particle(p, p->time_old, p->time, n, 1, 0);
		}
	      else
		output_a_particle(p, p->time_old, p->time, n, (int)(p->time == global_end_time), 0);
	      output_accum_ust += max(ust() - output_start_ust, 0.0);
	      */
	    }

	  p = p_tmp;
	}
    }

#ifdef PARALLEL
  if(Num_Proc > 0)
    {
      /*
	fprintf(stderr, "Proc %d wants to send %d particles.\n", ProcID, local_num_to_send);  fflush(stderr);
      */
      /* Now perform send/receives and advance_a_particle()'s for the
       * particles that are moving across processors. */
      MPI_Barrier(MPI_COMM_WORLD);
      num_to_send = 0;
      MPI_Allreduce(&local_num_to_send, &num_to_send, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      local_particle_transfers += local_num_to_send;
      if(!num_to_send)
	done = 1;
      else
	done = 0;

      zero_a_particle(&termination_p);
      termination_p.state = TERMINATION_MARKER;
      while(!done)
	{
	  /* This next loop does a full all processor to all processor
	   * exchange of particles.  This is done before any other moves are
	   * completed. */
	  for(i = 0; i < Num_Proc; i++)
	    {
	      if(i == ProcID) /* I'm the sender! */
		{
		  p = particles_to_send;
		  particles_to_send = NULL;
		  while(p)
		    {
		      p_tmp = p->next;
		      /*
			fprintf(stderr, "Proc %d is sending a particle to proc %d\n", i, p->owning_proc_id); fflush(stderr);
		      */
		      memcpy(&send_p, p, sizeof(particle_t));
		      mpi_retval = MPI_Send(&send_p, sizeof(particle_t), MPI_BYTE, p->owning_proc_id, 0, MPI_COMM_WORLD);
		      if(mpi_retval != MPI_SUCCESS)
			{
			  fprintf(stderr, "Proc %d failed on particle send (%d).\n", ProcID, mpi_retval);
			  exit(-1);
			}
		      free(p);
		      p = p_tmp;
		    }
		  for(j = 0; j < Num_Proc; j++)
		    if(i == j)
		      continue;
		    else
		      {
			/*
			  fprintf(stderr, "Proc %d is sending the TERMINATION_MARKER to proc %d\n", i, j); fflush(stderr);
			*/
			mpi_retval = MPI_Send(&termination_p, sizeof(particle_t), MPI_BYTE, j, 0, MPI_COMM_WORLD);
			if(mpi_retval != MPI_SUCCESS)
			  {
			    fprintf(stderr, "Proc %d failed on TERMINATION_MARKER send (%d).\n", ProcID, mpi_retval);
			    exit(-1);
			  }
		      }
		}
	      else
		{
		  recv_p.state = ACTIVE;
		  while(recv_p.state != TERMINATION_MARKER)
		    {
		      mpi_retval = MPI_Recv(&recv_p, sizeof(particle_t), MPI_BYTE, i, 0, MPI_COMM_WORLD, &mpi_status);
		      if(mpi_retval != MPI_SUCCESS)
			{
			  fprintf(stderr, "Proc %d failed on particle receive (%d).\n", ProcID, mpi_retval);
			  exit(-1);
			}
		      MPI_Get_count(&mpi_status, MPI_BYTE, &mpi_byte_count);
		  
		      /*
			fprintf(stderr, "Proc %d received a particle from proc %d with state = %d\n   MPI_SOURCE/TAG/ERROR = %d/%d/%d, MPI_Count = %d.\n", ProcID, i, recv_p.state, mpi_status.MPI_SOURCE, mpi_status.MPI_TAG, mpi_status.MPI_ERROR, mpi_byte_count);
			fflush(stderr);
		      */
		      if(mpi_byte_count != sizeof(particle_t))
			{
			  fprintf(stderr, "Proc %d did not have proper incoming byte count!\n", ProcID); fflush(stderr);
			  MPI_Finalize();
			  exit(-1);
			}

		      if(recv_p.state != TERMINATION_MARKER)
			{
			  add_to_do_list(&recv_p);
			  /*
			    fprintf(stderr, "Proc %d is adding a particle to its to do list.\n", ProcID); fflush(stderr);
			  */
			  /*
			    dump(NO_EXIT, &termination_p, "Rcvd by proc %d.\n", ProcID);
			  */
			}
		    }
		}
	      mpi_retval = MPI_Barrier(MPI_COMM_WORLD);
	      /*
		fprintf(stderr, "Proc %d BARRIER\n", ProcID);
	      */
	      if(mpi_retval != MPI_SUCCESS)
		{
		  fprintf(stderr, "Proc %d failed on MPI_BARRIER (%d).\n", ProcID, mpi_retval);
		  exit(-1);
		}
	    }

	  /* Now we move the new particles. */
	  local_num_to_send = 0;
	  p = particles_to_do;
	  particles_to_do = NULL;
	  while(p)
	    {
	      /*
		fprintf(stderr, "Proc %d is moving a particle.\n", ProcID); fflush(stderr);
	      */
	  
	      p_tmp = p->next;

	      /* This gets complicated as we need to simulate finishing
	       * one of the particle's little timesteps... */
	      if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
		dump1(EXIT, p, "New particle from other processor.");
	      /*
		dump(NO_EXIT, p, "Particle to be moved with new p->xi values.");
	      */
	  
	      /* Find the particles that have moved out of their element.
	       * Some of these have moved completely out of the computational
	       * domain. */
	      if(fabs(p->xi[0]) > 1.0 || fabs(p->xi[1]) > 1.0 || fabs(p->xi[2]) > 1.0)
		{
#ifdef SHOW_PARTICLE_MOVEMENT
		  fprintf(stderr, "PA");
#endif
		  find_exit_wound(p->owning_elem_id, p, p->x_old, p->x, p->xi, 0, 0);
		}

	      if(p->state == ACTIVE) /* <= 0 -> it left our processor or hit a wall. */
		{
		  /* Save the requested data for output later. */
		  if(p->output_sample_number != -1)
		    store_particle_data(p);
	  
		  /*
		    POSSIBLY NOT USED ANYMORE -- FROM CVS DIFF
		  output_start_ust = ust();
		  if(TimeIntegration == STEADY)
		    {
		      if(p->time == global_end_time)
			output_a_particle(p, p->old_time, p->time, n, 1, 0);
		    }
		  else
		    output_a_particle(p, p->old_time, p->time, n, (int)(p->time == global_end_time), 0);
		  output_accum_ust += MAX(ust() - output_start_ust, 0.0);
		  */
		  
		  advance_a_particle(p, global_start_time, global_end_time, n);
		}

	      if(p->state == PROC_TRANSFER) /* I need to go to yet another processor. */
		{
		  add_to_send_list(p);
		  local_num_to_send++;
		}
	      else if(p->state == DEAD) /* Delete me */
		free(p);	/* Bye bye */
	      else
		{
		  create_a_particle(p, p->owning_elem_id); /* This num_particles++ already */
		  free(p);
		}
	  
	      p = p_tmp;
	    }

	  num_to_send = 0;
	  MPI_Allreduce(&local_num_to_send, &num_to_send, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  local_particle_transfers += local_num_to_send;
	  if(!num_to_send)
	    done = 1;
	}
    }
#endif
  
  /* Call the routine that couples terms back to the continuum.  E.g.,
   * terms to be picked up in the momentum equation or a viscosity
   * value, etc... */ 
  couple_to_continuum();

  /* HERE , fix my timings, they are messed up somehow...  I might need
   * to send ust()'s between all the processors... An MPI_BARRIER
   * should do it, though. */
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  total_accum_ust = ust() - start_ust;
  particle_accum_ust = total_accum_ust - output_accum_ust;
#ifdef PARALLEL
  particle_accum_ust -= communication_accum_ust;
  local_num_particles = num_particles;
  num_particles = 0;
  particle_transfers = 0;
  local_max_particle_iterations = max_particle_iterations;
  /* OK to collect even if we ignore it... */
  local_max_newton_iterations = max_newton_iterations;
  local_particle_accum_ust = particle_accum_ust;
  local_output_accum_ust = output_accum_ust;
  temp_ust = ust();
  MPI_Reduce(&local_num_particles, &num_particles, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_max_particle_iterations, &max_particle_iterations, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_max_newton_iterations, &max_newton_iterations, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_particle_transfers, &particle_transfers, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_particle_accum_ust, &particle_accum_ust, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_output_accum_ust, &output_accum_ust, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  communication_accum_ust += MAX((dbl)ust() - (dbl)temp_ust, 0.0);
  local_communication_accum_ust = communication_accum_ust;
  MPI_Reduce(&local_communication_accum_ust, &communication_accum_ust, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  local_total_accum_ust = ust() - start_ust;
  MPI_Reduce(&local_total_accum_ust, &total_accum_ust, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
#endif
	  
  DPRINTF(stderr, "           Number of particles = %d\n", num_particles);
  DPRINTF(stderr, "       Max particle time steps = %d\n", max_particle_iterations);
  if(Particle_Model == INERTIAL_TRACER_IMPLICIT ||
     Particle_Model == SWIMMER_IMPLICIT ||
     Particle_Model == TRACER_IMPLICIT ||
     Particle_Model == CHARGED_TRACER_IMPLICIT ||
     Particle_Model == DIELECTROPHORETIC_TRACER_IMPLICIT)
    DPRINTF(stderr, "Max particle Newton iterations = %d\n", max_newton_iterations);
#ifdef PARALLEL
  DPRINTF(stderr, "            Particle transfers = %d\n", particle_transfers);
  DPRINTF(stderr, "               Particle ghosts = %d\n", total_num_particle_ghosts);
  DPRINTF(stderr, "            Communication time = %g seconds\n", communication_accum_ust);
#endif
  DPRINTF(stderr, "                 Particle time = %g seconds\n", particle_accum_ust);
  DPRINTF(stderr, "                   Output time = %g seconds\n", output_accum_ust);
  DPRINTF(stderr, "                    Total time = %g seconds\n", total_accum_ust);

  for(i = 0; i < Particle_Number_Sample_Types; i++)
    fflush(pa_fp[i]);
  if(Particle_Full_Output_Stride)
    fflush(pa_full_fp);

  return 0;
}


/* Fill the element_volume array with the volumes of each elements.
 * This can be called repeatedly (as in it doesn't allocate anything
 * internally).
 */
static dbl
fill_element_volumes(const int e_begin,
		     const int e_end)
{
  int i, ip, ip_total;
  dbl this_volume, total_volume, wt, det, xi[DIM];

  total_volume = 0.0;

  /* First get total element "volume" */
  for(i = e_begin; i < e_end; i++)
    {
#ifdef PARALLEL
      if(DPI_ptr->elem_owner[i] != ProcID)
	{
	  el_volume[i] = 0.0;
	  continue;
	}
#endif
      this_volume = 0.0;
      load_element_information(i);
      ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
      for(ip = 0; ip < ip_total; ip++)
	{
	  wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type);
	  if(pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
	     pd_glob[0]->CoordinateSystem == SWIRLING)
	    wt *= 2.0*M_PIE;

	  find_stu(ip, ei[pg->imtrx]->ielem_type, &(xi[0]), &(xi[1]), &(xi[2]));

	  /* Although I don't need to actually load_fv, and only need
	   * to load_basis_function() (I think), this routine is
	   * called once so I don't care about the cost.  Furthermore,
	   * I might need to do this (and more) when I put in
	   * deformable mesh. */
	  load_field_variables_at_xi(i, xi);

	  det = bf[pd->ShapeVar]->detJ * fv->h3;

	  this_volume += det * wt;
	}
      el_volume[i] = this_volume;
      total_volume += this_volume;
    }
  return total_volume;
}


/* Fills x with a uniformly distributed location in element el_index.
 * It does this in a possibly inefficient way, but it is only done
 * once, so who cares?  It uses a rejection sampling method.
 */
static int
position_particle_uniformly(const int elem_id,
			    particle_t * p)
{
  int i, j, done, counter;
  dbl min_x[DIM], max_x[DIM];

  memset(min_x, 0, DIM * sizeof(dbl));
  memset(max_x, 0, DIM * sizeof(dbl));
  load_element_node_coordinates(elem_id);
  for(i = 0; i < mdim; i++)
    {
      min_x[i] = node_coord[0][i];
      max_x[i] = node_coord[0][i];
    }
  for(i = 1; i < nodes_per_element; i++)
    for(j = 0; j < mdim; j++)
      {
	min_x[j] = MIN(min_x[j], node_coord[i][j]);
	max_x[j] = MAX(max_x[j], node_coord[i][j]);
      }
  
  done = 0;
  counter = 0;
  while(!done)
    {
      for(i = 0; i < mdim; i++)
	p->x[i] = min_x[i] + drand48() * (max_x[i] - min_x[i]);

      if(get_element_xi_newton(elem_id, p->x, p->xi) == -1)
	{
	  /* We cheat here and just call it a rejection if we cannot
	   * place the particle. */
#ifdef PARALLEL
	  fprintf(stderr, "Proc %d could not converge on element xi coordinates -- everything is OK\n", ProcID);
#else	
	  fprintf(stderr, "Could not converge on element xi coordinates -- everything is OK.\n");
#endif
	  return -1;
	}

      done = 1;
      for(i = 0; i < mdim; i++)
	if(fabs(p->xi[i]) >= 1.0)
	  done = 0;
      counter++;
      if(counter == 100 && !done)
	EH(-1, "Local element rejection sampling counter reached 100.");
    }
  return counter;
}


/* Just reminding myself that this will probably be very intersting
 * one day with output categorized by timestep or sample or
 * whatever...
 */
static char *
construct_filename(const char *template)
{
  static char result[256];

  strcpy(result, template);
#ifdef PARALLEL
  multiname(result, ProcID, Num_Proc);
#endif
  return result;
}


/* This routine solves for the location of intersection between a line
 * in physical coordinates and a boundary in element coordinates.  It
 * is a nonlinear Newton solve.
 *
 * Upon entrance xi_tmp contains the initial guess for Newton's method
 * in element coordinates, and fixed_index is coordinate index that is
 * fixed (and represents the edge/surface boundary).  The coeff array
 * is used in two different ways, dependent on a 2D or 3D usage.  In
 * 2D, coeff[0], coeff[1], and coeff[2] are the coefficients of the
 * line, coeff[0] * x + coeff[1] * y + coeff[2] = 0, and components 3,
 * 4, and 5 are ignored.  In 3D, the 0, 1, and 2 coefficients are
 * really the components of the "start" vector, x_start, and the 3, 4,
 * and 5 coefficents are the "direction" vector from x_start to x_end.
 * It's not much pre-calculation, but there you go...
 *
 * Upon exit, xi_tmp has the point of intersection.  It may have an
 * abs. value > 1.0, which indicates that the point of intersection
 * lies outside the element, and is probably not the intersection you
 * were looking for.
 * 
 * There is an obvious improvement to be done here, which is not to
 * fill all basis functions, but only the geometry interpolation
 * variables.  Similarly, don't go through load_fv, but just load the
 * spatial coordinates.  Also, if the dfdxi values are already
 * precomputed in the Jacobian, use them!
 *
 * Although it handles both 2D and 3D, it isn't done in a
 * "dimension-free" way.  This is used at a low level so I want it to
 * be fast. */
static void
get_boundary_xi_newton(const dbl * const coeff,
		       const int elem_id,
		       dbl * xi_tmp,
		       const int fixed_index,
		       const dbl error_tolerance)
{
  dbl update;
  dbl f_initial = 0.0, f, dfdxi1 = 0.0, t = 0.0, xi_norm = 0.0, f_scale = 0.0;
  dbl resid_initial = 0.0, resid_scale = 0.0;
  dbl pt[DIM], pt_xi[DIM];
  dbl tmp[DIM], tmp2[DIM];
  int i, si, iter, index1;
  struct Basis_Functions *map_bf;

  si = in_list(pd->IntegrationMap, 0, Num_Interpolations, Unique_Interpolations);
  map_bf = bfd[si];
  if(mdim == 2)
    {
      index1 = 1 - fixed_index;
    }
  else
    {
      index1 = (fixed_index + 1) % 3;
    }

  load_field_variables_at_xi(elem_id, xi_tmp);
  
  if(mdim == 2)
    {
      f = coeff[0] * fv->x[0] + coeff[1] * fv->x[1] + coeff[2];
      f_initial = fabs(f);
      f_scale = f_initial;
      xi_norm = fabs(xi_tmp[index1]);
      f_scale *= MAX(1.0, xi_norm);

      iter = 0;
      while(fabs(f) >= error_tolerance * f_scale &&
	    iter < 100)
	{
	  /*
	    fprintf(stderr, "  Candidate point is x = (%g,%g) and xi = (%g,%g), elem_id = %d, f = %g\n",
	    fv->x[0], fv->x[1], xi_tmp[0], xi_tmp[1], elem_id, f);
	  */

	  /* I think dfdxi might be available in another form, like
	   * bf[pd->ShapeVar]->J, but I need to check carefully.
	   * Checked (with Tom's help), and map_bf->J[a][b] =
	   * d_x[b]_d_xi[a].  Yay (it caused "no" speedup on the
	   * problem I happened to be playing with at the time, though
	   * it was not move-intensive...)
	   */
	  dfdxi1 = coeff[0] * map_bf->J[index1][0] + coeff[1] * map_bf->J[index1][1];
	  xi_tmp[index1] -= f/dfdxi1;
	  load_field_variables_at_xi(elem_id, xi_tmp);
	  
	  f = coeff[0] * fv->x[0] + coeff[1] * fv->x[1] + coeff[2];
	  f_scale = f_initial;
	  xi_norm = fabs(xi_tmp[index1]);
	  f_scale *= MAX(1.0, xi_norm);
	  
	  iter++;
	}
    }
  else
    {
      /*
      fprintf(stderr, "In _boundary_\n");
      */
      /*
      fprintf(stderr, "elem_id = %d\n", elem_id);
      fprintf(stderr, "coeff = (%g,%g,%g,%g,%g,%g)\n",
	      coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
      */
      memset(tmp, 0, DIM * sizeof(dbl));
      memcpy(tmp2, coeff, DIM * sizeof(dbl));
      if(get_element_xi_newton(elem_id, tmp2, tmp) == -1)
	EH(-1, "BAD");
      t = 1.0;
      for(i = 0; i < DIM; i++)
	pt[i] = coeff[i] + t * coeff[i+3];
      memset(pt_xi, 0, DIM * sizeof(dbl));
      if(get_element_xi_newton(elem_id, pt, pt_xi) == -1)
	{
	  fprintf(stderr, "Bad gexn(), pt = (%.10g,%.10g,%.10g), pt_xi = (%.10g,%.10g,%.10g)\n",
		  pt[0], pt[1], pt[2], pt_xi[0], pt_xi[1], pt_xi[2]);
	}
      /*
      fprintf(stderr, "tmp(0.0) = (%.10g,%.10g,%.10g), pt_xi(1.0) = (%.10g,%.10g,%.10g)\n",
	      tmp[0], tmp[1], tmp[2], pt_xi[0], pt_xi[1], pt_xi[2]);
      fprintf(stderr, "t = (%.10g-%.10g)/(%.10g-%.10g) = %g\n",
	      xi_tmp[fixed_index], tmp[fixed_index], pt_xi[fixed_index], tmp[fixed_index],
	      (xi_tmp[fixed_index] - tmp[fixed_index]) / (pt_xi[fixed_index] - tmp[fixed_index]));
      */
      t = (xi_tmp[fixed_index] - tmp[fixed_index]) / (pt_xi[fixed_index] - tmp[fixed_index]);
      for(i = 0; i < DIM; i++)
	pt[i] = coeff[i] + t * coeff[i+3];
      memset(pt_xi, 0, DIM * sizeof(dbl));
      if(get_element_xi_newton(elem_id, pt, pt_xi) == -1)
	EH(-1, "ALSO BAD");
      /*
      fprintf(stderr, "t = %g, pt = (%.10g,%.10g,%.10g), pt_xi = (%.10g,%.10g,%.10g)\n",
	      t, pt[0], pt[1], pt[2], pt_xi[0], pt_xi[1], pt_xi[2]);
      */
      f = pt_xi[fixed_index] - xi_tmp[fixed_index];
      resid_scale = resid_initial = MAX(fabs(f), 1.0);
      xi_norm = nnorm(3, pt_xi);
      resid_scale *= MAX(1.0, xi_norm);

      /*
      fprintf(stderr, "STARTING\n");
      fprintf(stderr, "fv->x = (%.10g,%.10g,%.10g)\npt = (%.10g,%.10g,%.10g)\npt_xi = (%.10g,%.10g,%.10g)\ncoeff = (%.10g,%.10g,%.10g,%.10g,%.10g,%.10g)\nt = %.10g\nresid = %.10g\n",
	      fv->x[0], fv->x[1], fv->x[2],
	      pt[0], pt[1], pt[2],
	      pt_xi[0], pt_xi[1], pt_xi[2],
	      coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5],
	      t, f);
      fprintf(stderr, "Termination value = %.10g\n", 1.0e-8 * resid_scale);
      */
      
      iter = 0;
      while(fabs(f) >= error_tolerance * resid_scale &&
	    iter < 100)
	{
	  /*
	  fprintf(stderr, "iter = %d\n", iter);
	  fprintf(stderr, "  starting t = %.10g, starting f = %.10g\n", t, f);
	  */
	  dfdxi1 = 0.0;
	  for(i = 0; i < DIM; i++)
	    dfdxi1 += map_bf->B[i][fixed_index] * coeff[i+3];
	  update = f/dfdxi1;
	  t -= update;

	  for(i = 0; i < DIM; i++)
	    pt[i] = coeff[i] + t * coeff[i+3];
	  memset(pt_xi, 0, DIM * sizeof(dbl));
	  if(get_element_xi_newton(elem_id, pt, pt_xi) == -1)
	    {
	      fprintf(stderr, "  dfdxi1 = %.10g, update = %.10g\n",
		      dfdxi1, update);
	      fprintf(stderr, "  pt = (%.10g,%.10g,%.10g), pt_xi = (%.10g,%.10g,%.10g)\n",
		      pt[0], pt[1], pt[2], pt_xi[0], pt_xi[1], pt_xi[2]);
	      fprintf(stderr, "  resid_scale = %.10g, resid_initial = %.10g\n",
		      resid_scale, resid_initial);
	      fprintf(stderr, "  coeff = (%.10g,%.10g,%.10g,%.10g,%.10g,%.10g)\n",
		      coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
	    }
	  f = pt_xi[fixed_index] - xi_tmp[fixed_index];
	  xi_norm = nnorm(3, pt_xi);
	  resid_scale = resid_initial;
	  resid_scale *= MAX(1.0, xi_norm);

	  /*
	  fprintf(stderr, "  dfdxi1 = %.10g, update = %.10g\n",
		  dfdxi1, update);
	  fprintf(stderr, "  pt = (%.10g,%.10g,%.10g), pt_xi = (%.10g,%.10g,%.10g)\n",
		  pt[0], pt[1], pt[2], pt_xi[0], pt_xi[1], pt_xi[2]);
	  fprintf(stderr, "  resid_scale = %.10g, resid_initial = %.10g\n",
		  resid_scale, resid_initial);
	  fprintf(stderr, "  coeff = (%.10g,%.10g,%.10g,%.10g,%.10g,%.10g)\n",
		  coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
	  */
	  
	  iter++;
	}
      memcpy(xi_tmp, pt_xi, DIM * sizeof(dbl));
    }

  if(iter == 100 && fabs(f) > error_tolerance)
    {
      fprintf(stderr, "last dfdxi = %g, t = %g, f = %g, resid_scale = %g, resid_initial = %g, xi_norm = %g, pt = (%g,%g,%g), pt_xi = (%g,%g,%g)\n",
	      dfdxi1, t, f, resid_scale, resid_initial, xi_norm, pt[0], pt[1], pt[2], pt_xi[0], pt_xi[1], pt_xi[2]);
      EH(-1, "Exceeded 100 Newton iterations in get_boundary_xi_newton()");
    }
     
  /*
  fprintf(stderr, "Exiting at x = (%g,%g,%g), xi = (%g,%g,%g), resid/f = %g\n",
	 fv->x[0], fv->x[1], fv->x[2], xi_tmp[0], xi_tmp[1], xi_tmp[2], (mdim==2?f:resid_norm));
  */
}


/* Create a new particle_t item in the current particle_head headed
 * link list.  This routine needs an elem_id to properly insert the
 * new particle.  It sets the p->owning_elem_id to this in case you
 * forgot. It returns a pointer to the new particle, with the next and
 * last fields filled.
 */
static particle_t *
obtain_particle_space(const int elem_id)
{
  particle_t * p;
  char s[128];
  
  if(!element_particle_list_head[elem_id])
    {
      element_particle_list_head[elem_id] = (particle_t *)malloc(sizeof(particle_t));
      if(!element_particle_list_head[elem_id])
	{
	  sprintf(s, "Could not allocate %ld bytes of space for particle.\n", sizeof(particle_t));
	  EH(-1, s);
	}
      p = element_particle_list_head[elem_id];
      p->last = NULL;
      p->next = NULL;
    }
  else
    {
      p = (particle_t *)malloc(sizeof(particle_t));
      if(!p)
	EH(-1, "Could not allocate space for particle.");
      p->next = element_particle_list_head[elem_id];
      p->last = NULL;
      element_particle_list_head[elem_id]->last = p;
      element_particle_list_head[elem_id] = p;
    }
  p->owning_elem_id = elem_id;
  
  return p;
}


/* This routine creates a particle.  It seeds this new particle with
 * the passed one.  It returns the new particle.
 *
 * It also enforces any Particle_Move_Domain data.  It should also enforce
 * PBC_IMPERMEABLE's, too...
 */
static particle_t *
create_a_particle(particle_t * seed_p,
		  const int elem_id)
{
  particle_t * p, *p_last_keep, *p_next_keep;
  int i;

  if(DPI_ptr->elem_owner[elem_id] != ProcID &&
     seed_p->state != GHOST)
    dump4(EXIT, seed_p, "Proc %d tried to create a particle on element %d(%d), which it doesn't own.\n", ProcID, elem_id, DPI_ptr->elem_index_global[elem_id]);
  
  p = obtain_particle_space(elem_id);
  p_last_keep = p->last;
  p_next_keep = p->next;
  memcpy(p, seed_p, sizeof(particle_t));
  p->last = p_last_keep;
  p->next = p_next_keep;
  
  /* Enforce particle bounds.  Heavy-handed, but effective...  This
   * should be "OK" because the bounds should be "very close" to the
   * actual boundaries so we're cutting off a "negligable" amount of
   * volume.  In fact, this shouldn't be hit anymore since I modified
   * the rejection sampling method to handle this.  Leaving it here
   * "just in case". */
  if(Particle_Move_Domain == BRICK)
    for(i = 0; i < pdim; i++)
      {
	p->x[i] = MAX(p->x[i], Particle_Move_Domain_Reals[2*i]);
	p->x[i] = MIN(p->x[i], Particle_Move_Domain_Reals[2*i+1]);
      }
  
  num_particles++;
  return p;
}


/* This routine performs a rejection sampling method to create a
 * particle at the initialization step.  It is possible to specify
 * particle creation bounds in such a way that *all* particles are
 * rejected, to be careful.  An "always rejected" processor domain
 * problem in parallel is avoided by controlling the creation of
 * particles at the globla (processor 0) level.
 *
 * There are two ways the particle can be rejected.  If the local
 * element is non-rectilinear and the local sampling creates a
 * particle outside the element (fabs(xi[i])>1.0), then the particle
 * is rejected, even though its position might lie in an adjacent
 * element just fine.  This occurs in position_particle_uniformly, and
 * we don't need to worry about it here.  The particle may also be
 * rejected on a global basis depending on it's location and the
 * enforced particle bounds.
 */
static int
rejection_sample_a_particle(void)
{
  particle_t p;
  int i, el_index, rejection, num_element_samples;
  dbl r, current_volume;

  r = drand48() * my_volume;

  /* This method of selecting an element to introduce a particle is
   * very sensitive to the DPI_ptr->elem_owner[i] array.  el_volume[i]
   * = 0 if that element is not owned by this processor, so don't try
   * to create a particle there! */
  el_index = 0;
  current_volume = 0.0;
  while(current_volume < r)
    current_volume += el_volume[el_index++];
  el_index--;
  
  /* return value is currently ignored... */
  zero_a_particle(&p);
  num_element_samples = position_particle_uniformly(el_index, &p);
  if(num_element_samples == -1)
    rejection = 1;
  else
    {
      rejection = 0;
      if(Particle_Creation_Domain == BRICK)
	{
	  for(i = 0; i < pdim; i++)
	    if(p.x[i] < Particle_Creation_Domain_Reals[2*i] || p.x[i] > Particle_Creation_Domain_Reals[2*i+1])
	      rejection = 1;
	}
    }

  if(!rejection)
    initialize_and_create_new_particle(&p, el_index);

  return rejection;
}


/* Assuming you just created a new particle (not ghosted, and
 * generally not from the boundary inflow condition), this is the
 * routine that will initialize the particle's data (e.g., velocity),
 * and call the create_a_particle function to place it in the element
 * lists.  The particle's p.x value must be filled in and correct.  */
static void
initialize_and_create_new_particle(particle_t *p,
				   const int elem_id)
{
  memcpy(p->x_old, p->x, DIM * sizeof(dbl));
  memset(p->xi, 0, DIM * sizeof(dbl));
  if(get_element_xi_newton(elem_id, p->x, p->xi) == -1)
    dump1(EXIT, p, "Could not find element coordinates.");
  p->owning_elem_id = elem_id;
  p = create_a_particle(p, elem_id);
  load_field_variables_at_xi(p->owning_elem_id, p->xi);
  if(TimeIntegration == STEADY)
    p->time = 0.0;
  else
    p->time = tran->init_time;
  
  switch(Particle_Model)
    {
    case INERTIAL_TRACER_EXPLICIT:
    case INERTIAL_TRACER_IMPLICIT:
    case TRACER_EXPLICIT:
    case TRACER_IMPLICIT:
    case CHARGED_TRACER_EXPLICIT:
    case CHARGED_TRACER_IMPLICIT:
    case DIELECTROPHORETIC_TRACER_IMPLICIT:
      memcpy(p->v, fv->v, DIM * sizeof(dbl));
      break;
    case SWIMMER_IMPLICIT:
    case SWIMMER_EXPLICIT:
      p->theta = p->theta_old = drand48()*2.0*M_PIE; /* theta=pi => pointing up! */
      if(pdim == 3)
	p->phi = p->phi_old = 0.0;
      memcpy(p->v, fv->v, DIM * sizeof(dbl));
      break;
    default:
      dump1(EXIT, p, "Unknown Particle_Model.  You shouldn't be here.");
    }
  if(p->output_sample_number > -1)
    store_particle_data(p);
  output_a_particle(p, 0.0, 0.0, 0, 1, 1);
}


/* This routine will create a (possible large) structure that reverse
 * maps element sides to PBC's in the input.  We need to do this for
 * fast particle <-> surface interaction, sorry.
 *
 * This now includes a value for an off-processor element.  The value
 * is -2 + owning processor ID (results <= -2).  -1 still indicates an
 * actual domain boundary.
 */
static void
create_element_particle_info_maps(void)
{
  int i, j;
  int SS_id, PBC_id, elem_id, side_id;
#ifdef PARALLEL
  int proc_id, neighbor_elem_id, global_elem_id;
  int node_id, local_node_id, global_node_id, ghost_id;
  int indices[4], node_compares[4], elem_ptr_start, nodes_per_side, match_found;
  int msg_tag, msg_size, *msg, tmp;
  MPI_Status mpi_status;

  int output = 0;
  char fp2name[80];
  FILE *fp, *fp2;
  /*  static void test_map_integrity(void); */

  /* initialize indices */
  for (i = 0; i < 4; i++) {
    indices[i] = 0;
  }
#endif

  /* Assuming that you will not apply more than one PBC to a sideset! */
  for(i = 0; i < static_exo->num_side_sets; i++)
    {
      /* Get this sidesets name (id) */
      SS_id = static_exo->ss_id[i];

      /* Is that one we're using?  If so, save the PBC index
       * corresponding to it. */
      PBC_id = -1;
      for(j = 0; j < Particle_Number_PBCs && PBC_id == -1; j++)
	if(PBCs[j].SS_id == SS_id)
	  PBC_id = j;
      if(PBC_id == -1)
	continue;

      for(j = 0; j < static_exo->ss_num_sides[i]; j++)
	{
	  elem_id = static_exo->ss_elem_list[static_exo->ss_elem_index[i] + j];
	  side_id = static_exo->ss_side_list[static_exo->ss_elem_index[i] + j];
	  /* To get 0-based side id's */
	  side_id--;
	  /*
	  fprintf(stderr, "Setting element %d side %d to %d\n",
		 elem_id, side_id, PBC_id);
	  */
	  element_particle_info[elem_id].PBC_side_id[side_id] = PBC_id;
	}
    }

#ifdef PARALLEL
  /* This part gets tricky.  For the boundary elements, there are
   * possibly > 1 processors that have a copy of them.  We need to
   * indicate to everyone who has copies and correlate the elem <->
   * elem mapping.  Furthermore, ghost particles need to exist on the
   * other processors so the proper nodal contributions due to
   * particles are calculated.  Great.
   *
   * Outgoing broadcasts:
   * -1 global_elem_id local_elem_id => I'm a ghost and need to tell someone.
   * global_elem_id global_node_id_1 ... global_node_id_N I need a
   *              matching element on the other side from another proc.
   * -2 => I'm done, go to next proc.
   */
  if(Num_Proc == 1)
    return; /* Nothing to do, right? */

  nodes_per_side = (mdim == 2 ? 2 : 4);
  msg_size = nodes_per_side + 1; /* num side nodes + 1 */
  msg = (int *)calloc((unsigned)msg_size, sizeof(int));
  msg_tag = 0;
  if(output)
    {
      sprintf(fp2name, "map%d_response.txt", ProcID);
      fp2 = fopen(fp2name, "w");
    }
  for(proc_id = 0; proc_id < Num_Proc; proc_id++)
    {
      DPRINTF(stderr, "%d, ", proc_id); fflush(stderr);
      if(proc_id == ProcID)
	{
	  if(output)
	    {
	      fp = fopen("map_creation.txt", proc_id?"a":"w");
	      if(!fp)
		EH(-1, "Could not open map_creation.txt");
	    }

	  for(elem_id = 0; elem_id < static_exo->num_elems; elem_id++)
	    {	    
	      global_elem_id = DPI_ptr->elem_index_global[elem_id];

	      /* If we're a ghosted element, broadcast so someone
	       * knows they should send us ghost information. */
	      if(DPI_ptr->elem_owner[elem_id] != proc_id)
		{
		  msg[0] = -1;
		  msg[1] = global_elem_id;
		  msg[2] = elem_id;
		  if(output)
		    {
		      fprintf(fp, "Proc %d Elem %4d(%3d) not owned, sending ghost message.\n", ProcID, global_elem_id, elem_id);
		      fflush(fp);
		    }
		  MPI_Barrier(MPI_COMM_WORLD);
		  MPI_Bcast(msg, msg_size, MPI_INT, proc_id, MPI_COMM_WORLD);
		}
	      /* Otherwise we're looking for matched global node ids
	       * on unmatched surfaces. */
	      else
		for(side_id = 0; side_id < sides_per_element; side_id++)
		  {
		    neighbor_elem_id = static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + side_id];
		    if(element_particle_info[elem_id].PBC_side_id[side_id] == -1 &&
		       (neighbor_elem_id == -1 || DPI_ptr->elem_owner[neighbor_elem_id] != ProcID))
		      {
			msg[0] = global_elem_id;
			if(mdim == 2)
			  {
			    indices[0] = side_id;
			    indices[1] = (side_id + 1) % 4;
			  }
			else
			  fill_hex_side_indices(side_id, indices);
			
			/* I want global node id's, not node numbers here! */
			elem_ptr_start = Proc_Connect_Ptr[elem_id];
			for(node_id = 0; node_id < nodes_per_side; node_id++)
			  {
			    local_node_id = Proc_Elem_Connect[elem_ptr_start + indices[node_id]];
			    global_node_id = DPI_ptr->node_index_global[local_node_id];
			    msg[node_id+1] = global_node_id;
			  }
			
			/* I thought about sorting the node ids for
			 * some sort of faster calculation but I
			 * didn't think there was nothing to gain as
			 * long as sender/receiver generate the list
			 * in the exact same way.  Now I think I need
			 * to sort them because elements on different
			 * procs might order them differently... */
			/* Bubble sort.  Ugh. */
			for(i = 0; i < nodes_per_side - 1; i++)
			  for(j = i + 1; j < nodes_per_side; j++)
			    if(msg[i+1] > msg[j+1])
			      {
				tmp = msg[i+1];
				msg[i+1] = msg[j+1];
				msg[j+1] = tmp;
			      }

			if(output)
			  {
			    fprintf(fp, "Proc %d Elem %4d(%3d) Side %d not linked nodes =", ProcID, global_elem_id, elem_id, side_id);
			    for(i = 0; i < nodes_per_side; i++)
			      fprintf(fp, " %5d", msg[i+1]);
			    fprintf(fp, "\n");
			    fflush(fp);
			  }


			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(msg, msg_size, MPI_INT, proc_id, MPI_COMM_WORLD);
			MPI_Recv(msg, 2, MPI_INT, MPI_ANY_SOURCE, msg_tag, MPI_COMM_WORLD, &mpi_status);
			if(output)
			  {
			    fprintf(fp, "\tproc %d claimed the other side as its local element %3d\n", msg[0], msg[1]);
			    fflush(fp);
			  }
			element_particle_info[elem_id].PBC_side_id[side_id] = -(2 + msg[0]);
			element_particle_info[elem_id].owner_local_element_id[side_id] = msg[1];
			msg_tag++;
		      } /* if unknown reference ... */
		  } /* for side_id ... */
	    } /* for elem_id ... */
	  if(output)
	    {
	      fprintf(fp, "Proc %d sending termination ticket.\n", ProcID);
	      fclose(fp);
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	  msg[0] = -2;		/* Send termination */
	  MPI_Bcast(msg, msg_size, MPI_INT, proc_id, MPI_COMM_WORLD);
	} /* if proc_id == ProcID */
      else
	{
	  msg[0] = 0;
	  while(msg[0] != -2)
	    {
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Bcast(msg, msg_size, MPI_INT, proc_id, MPI_COMM_WORLD);
	      if(output)
		{
		  fprintf(fp2, "Received ");
		  for(i = 0; i < msg_size; i++)
		    fprintf(fp2, " %5d", msg[i]);
		  fprintf(fp2, "\n");
		  fflush(fp2);
		}
	      if(msg[0] == -2)
		{
		  if(output)
		    {
		      fprintf(fp2, "\tTerminating loop for proc %d\n", proc_id);
		      fflush(fp2);
		    }
		  break;
		}
	      if(msg[0] == -1)
		{
		  global_elem_id = msg[1];
		  if(output)
		    {
		      fprintf(fp2, "\tLooking for ghosting elem %4d ... ", global_elem_id);
		      fflush(fp2);
		    }
		  elem_id = in_list(global_elem_id, 0, static_exo->num_elems, DPI_ptr->elem_index_global);
		  if(elem_id > -1 && DPI_ptr->elem_owner[elem_id] == ProcID)
		    {
		      ghost_id = element_particle_info[elem_id].num_ghost_target_elems++;
		      element_particle_info[elem_id].ghost_proc = (int *)realloc(element_particle_info[elem_id].ghost_proc, (ghost_id + 1) * sizeof(int));
		      element_particle_info[elem_id].ghost_local_elem_id = (int *)realloc(element_particle_info[elem_id].ghost_local_elem_id, (ghost_id + 1) * sizeof(int));
		      element_particle_info[elem_id].ghost_proc[ghost_id] = proc_id;
		      element_particle_info[elem_id].ghost_local_elem_id[ghost_id] = msg[2];
		      if(output)
			{
			  fprintf(fp2, "from local elem %3d to proc=%d local_elem=%3d\n", elem_id, proc_id, msg[2]);
			  fflush(fp2);
			}
		    }
		  else
		    {
		      if(output)
			{
			  fprintf(fp2, "don't own it.\n");
			  fflush(fp2);
			}
		    }
		}
	      else
		{
		  global_elem_id = msg[0];
		  if(output && 0)
		    {
		      fprintf(fp2, "\tSearching for nodes");
		      for(i = 0; i < nodes_per_side; i++)
			fprintf(fp2, " %5d", msg[i+1]);
		      fprintf(fp2, " ... ");
		      fflush(fp2);
		    }
		  match_found = 0;
		  for(elem_id = 0; elem_id < static_exo->num_elems && !match_found; elem_id++)
		    {
		      if(DPI_ptr->elem_owner[elem_id] != ProcID)
			continue;
		      elem_ptr_start = Proc_Connect_Ptr[elem_id];
		      for(side_id = 0; side_id < sides_per_element && !match_found; side_id++)
			{
			  /* If the element on the other side is owned by us, skip this side. */
			  neighbor_elem_id = static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + side_id];
			  if(DPI_ptr->elem_owner[neighbor_elem_id] == ProcID)
			    continue;
			  
			  if(mdim == 2)
			    {
			      indices[0] = side_id;
			      indices[1] = (side_id + 1) % 4;
			    }
			  else
			    fill_hex_side_indices(side_id, indices);

			  for(node_id = 0; node_id < nodes_per_side; node_id++)
			    {
			      local_node_id = Proc_Elem_Connect[elem_ptr_start + indices[node_id]];
			      global_node_id = DPI_ptr->node_index_global[local_node_id];
			      node_compares[node_id] = global_node_id;
			    }

			  /* Bubble sort.  Ugh. */
			  for(i = 0; i < nodes_per_side - 1; i++)
			    for(j = i + 1; j < nodes_per_side; j++)
			      if(node_compares[i] > node_compares[j])
				{
				  tmp = node_compares[i];
				  node_compares[i] = node_compares[j];
				  node_compares[j] = tmp;
				}
			  
			  if(output && 0)
			    {
			      fprintf(fp2, "\tComparing to %d %d %d %d\n", node_compares[0], node_compares[1], node_compares[2], node_compares[3]);
			      fflush(fp2);
			    }
			  
			  /* If they share *any* nodes */
			  /*
			  match_found = 0;
			  for(node_id = 0; node_id < nodes_per_side && !match_found; node_id++)
			    for(i = 0; i < nodes_per_side && !match_found; i++)
			      if(node_compares[node_id] == msg[i+1])
			      match_found = 1;
			  */
			  
			  /* If they share a face */

			  match_found = 1;
			  for(node_id = 0; node_id < nodes_per_side && match_found; node_id++)
			    if(node_compares[node_id] != msg[node_id + 1])
			      match_found = 0;
			} /* for side_id ... */
		    } /* for elem_id ... */
		  if(match_found)
		    {
		      elem_id--;
		      if(output)
			{
			  fprintf(fp2, "matching elem_id = %d\n", elem_id);
			  fflush(fp2);
			  fprintf(fp2, " found at elem %4d(%3d)\n", DPI_ptr->elem_index_global[elem_id], elem_id);
			  fflush(fp2);
			}
		      msg[0] = ProcID;
		      msg[1] = elem_id;
		      MPI_Send(msg, 2, MPI_INT, proc_id, msg_tag, MPI_COMM_WORLD);
		    }
		  else
		    {
		      if(output)
			{
			  fprintf(fp2, " not found.\n");
			  fflush(fp2);
			}
		    }
		  msg_tag++;
		} /* msg[0] indicates node matching */
	    } /* while(not terminated) */
	} /* proc_id != ProcID */
    } /* for i = 0 ... Num_Proc-1 */
  free(msg);
  if(output)
    fclose(fp2);

  test_map_integrity();
#endif
}


/* This routine determines what happens when a particle intersects a
 * surface on its way out of the domain.  The intersection was already
 * determined earlier. */
static void
handle_surface_interaction(particle_t * p,
			   dbl * x,
			   dbl * xi,
			   int PBC_id)
{
  int i;
  PBC_t *PBC;
  FILE *fp;

  PBC = &(PBCs[PBC_id]);

  switch(PBC->type)
    {
    case PBC_IMPERMEABLE:
      fprintf(stderr, "PARTICLE ACTUALLY REACHED IMPERMEABLE BOUNDARY\n");
      /* fall through */ /* And follow through to be deleted "normally". */
      /* fall through */
    case PBC_OUTFLOW:
    case PBC_SOURCE:		/* this is not quite right, but oh well... */
    case PBC_FREESTREAM_SOURCE:	/* this is not quite right, but oh well... */
      /*
      fprintf(stderr, "Proc %d killed particle (1)\n", ProcID); fflush(stderr);
      */
      p->state = DEAD;
      break;
    case PBC_TARGET:
      /*
      fprintf(stderr, "Proc %d killed particle (2)\n", ProcID); fflush(stderr);
      */
      fp = fopen(construct_filename(PBC->string_data), "a");
      for(i = 0; i < pdim; i++)
	fprintf(fp, "%g ", x[i]);
      fprintf(fp, "\n");
      fclose(fp);
      p->state = DEAD;
      break;
    default:
      dump3(EXIT, p, "Unknown PBC type (%d) found on PBC_id %d.\n", PBC->type, PBC_id);
    }
}


/* This routine introduces new particles due to SOURCE particle
 * boundary conditions (PBC's).  It does so assuming the source SS's
 * side is linear.  Shouldn't be too big a deal to assume that, but
 * one day I need to fix it (for quadratic deformable mesh and
 * CYLINDRICAL, too).
 */
static void
generate_source_particles(const dbl tt,	/* parameter to vary time integration */
			  const dbl dt, /* time step size */
			  const dbl creation_time) /* time when particle was created. */
{
  PBC_t *PBC;
  particle_t p;
  int SS_id, elem_id, side_id, ip_id;
  int i, j, k, m, ip_total, num_nodes_on_side;
  int *local_ss_node_list, local_elem_node_id[MAX_NODES_PER_SIDE];
  dbl inflow_volume, expected_particles, particle_volume, xi[DIM], s, t, u, wt;

  for(i = 0; i < Particle_Number_PBCs; i++)
    {
      PBC = &(PBCs[i]);
      if(PBC->type == PBC_SOURCE ||
	 PBC->type == PBC_FREESTREAM_SOURCE)
	{
	  SS_id = PBC->SS_id;
	  for(j = 0; j < static_exo->num_side_sets; j++)
	    if(SS_id == static_exo->ss_id[j])
	      {
		for(k = 0; k < static_exo->ss_num_sides[j]; k++)
		  {

		    elem_id = static_exo->ss_elem_list[static_exo->ss_elem_index[j] + k];
		    side_id = static_exo->ss_side_list[static_exo->ss_elem_index[j] + k];
		    side_id--;	/* make it 0-based */

#ifdef PARALLEL
		    /* Don't try to create particles for elements that
		     * don't belong on this processor!. */
		    if(DPI_ptr->elem_owner[elem_id] != ProcID)
		      continue;
#endif
		    
		    if(pd->CoordinateSystem != CARTESIAN &&
		       pd->CoordinateSystem != CYLINDRICAL &&
		       pd->CoordinateSystem != SWIRLING)
		      EH(-1, "Cannot apply volumetric source particles on this coordinate system.");

		    load_element_information(elem_id);
		    num_nodes_on_side = static_exo->ss_node_side_index[j][k+1] - static_exo->ss_node_side_index[j][k];
		    local_ss_node_list = &(static_exo->ss_node_list[j][static_exo->ss_node_side_index[j][k]]);
		    for(m = 0; m < num_nodes_on_side; m++)
		      {
			if((local_elem_node_id[m] = in_list(local_ss_node_list[m], 0, nodes_per_element, &(static_exo->elem_node_list[ei[pg->imtrx]->iconnect_ptr]))) == -1)
			  EH(-1, "Bad id_local_elem_coord lookup");
		      }

		    /* For Q1 elements, the node list I get from above is the same as:
		       fill_hex_side_indices(side_id, tmp_indices);
		    * I expect it would be different for Q2...
		    */
		    
		    /* First we want to get the inflow volume. */
		    inflow_volume = 0.0;
		    ip_total = elem_info(NQUAD_SURF, ei[pg->imtrx]->ielem_type);
		    for(ip_id = 0; ip_id < ip_total; ip_id++)
		      {
			/* I don't think I want it zero-based according to the comment in rf_bc_const.h */
			find_surf_st(ip_id, ei[pg->imtrx]->ielem_type, side_id + 1, pd->Num_Dim, xi, &s, &t, &u);
			wt = Gq_surf_weight(ip_id, ei[pg->imtrx]->ielem_type);
			if(pd->CoordinateSystem == CYLINDRICAL ||
			   pd->CoordinateSystem == SWIRLING)
			  wt *= 2.0*M_PIE;
			load_field_variables_at_xi(elem_id, xi);
			surface_determinant_and_normal(elem_id, ei[pg->imtrx]->iconnect_ptr, ei[pg->imtrx]->num_local_nodes, 
						       ei[pg->imtrx]->ielem_dim - 1, side_id + 1,
						       num_nodes_on_side, local_elem_node_id);
			for(m = 0; m < mdim; m++)
			  inflow_volume += -(fv->v[m] * fv->snormal[m]) * wt * fv->sdet; /* Normal is outward pointing. */
		      }
		    inflow_volume *= dt;
		    particle_volume = 4.0/3.0 * M_PIE * pow(Particle_Radius, 3.0);
		    expected_particles = inflow_volume * PBC->real_data[0] / particle_volume;
		    
		    if(expected_particles < 0.0)
		      {
			fprintf(stderr, "# expected particles = %g on elem %d side %d with inflow volume = %g\n",
				expected_particles, elem_id, side_id, inflow_volume);
			EH(-1, "Negative number of expected particles.");
		      }

		    for(m = 0; m < floor(expected_particles); m++)
		      {
			zero_a_particle(&p);
			select_random_side_location(elem_id, side_id, p.x, p.xi);
			memcpy(p.xi_old, p.xi, DIM*sizeof(dbl));
			memcpy(p.x_old, p.x, DIM*sizeof(dbl));
			p.owning_elem_id = elem_id;
			if(PBC->type == PBC_FREESTREAM_SOURCE)
			  {
			    load_field_variables_at_xi(elem_id, p.xi);
			    memcpy(p.v, fv->v, DIM * sizeof(dbl));
			    memcpy(p.v_old, p.v, DIM * sizeof(dbl));
			  }
			p.time = creation_time;
			create_a_particle(&p, elem_id);
		      }
		    if(drand48() < (expected_particles - floor(expected_particles)))
		      {
			zero_a_particle(&p);
			select_random_side_location(elem_id, side_id, p.x, p.xi);
			memcpy(p.xi_old, p.xi, DIM*sizeof(dbl));
			memcpy(p.x_old, p.x, DIM*sizeof(dbl));
			p.owning_elem_id = elem_id;
			if(PBC->type == PBC_FREESTREAM_SOURCE)
			  {
			    load_field_variables_at_xi(elem_id, p.xi);
			    memcpy(p.v, fv->v, DIM * sizeof(dbl));
			    memcpy(p.v_old, p.v, DIM * sizeof(dbl));
			  }
			p.time = creation_time;
			create_a_particle(&p, elem_id);
		      }
		  } /* for k ... static_exo->ss_num_sides[j] */
	      } /* SS_id == static_exo->ss_id[j] */
	} /* PBC_Type == PBC_SOURCE || PBC_FREESTREAM_SOURCE */
    } /* i = ... Particle_Number_PBCs */
}


/* This routine will initialize the parameters in a sane way to the
 * passed particle.*/
static void
zero_a_particle(particle_t *p)
{
  p->state = ACTIVE;
  p->theta = M_PIE;
  p->theta_old = M_PIE;
  p->output_sample_number = -1;
  p->owning_elem_id = -1;
  p->time = 0.0;
  p->time_old = 0.0;
  p->next = NULL;
  p->last = NULL;
  memset(p->real_data, 0, MAX_DATA_REAL_VALUES * sizeof(dbl));
  memset(p->real_data_old, 0, MAX_DATA_REAL_VALUES * sizeof(dbl));
  memset(p->int_data, 0, MAX_DATA_INT_VALUES * sizeof(int));
  memset(p->int_data_old, 0, MAX_DATA_INT_VALUES * sizeof(int));
  memset(p->x, 0, DIM * sizeof(dbl));
  memset(p->x_old, 0, DIM * sizeof(dbl));
  memset(p->xi, 0, DIM * sizeof(dbl));
  memset(p->xi_old, 0, DIM * sizeof(dbl));
  memset(p->v, 0, DIM * sizeof(dbl));
  memset(p->v_old, 0, DIM * sizeof(dbl));
#ifdef PARALLEL
  memset(p->x_start, 0, DIM * sizeof(dbl));
  memset(p->x_end, 0, DIM * sizeof(dbl));
  p->owning_proc_id = ProcID;
#endif

}



/* This routine is called from the main particle initialization
 * routine.  It's job currently is just to initialize a file to log
 * particle intersection locations. */
static void
initialize_surface_interactions(void)
{
  PBC_t *PBC;
  FILE *fp;
  int i;
  char err_msg[128];

  Num_Impermeable_PBCs = 0;
  for(i = 0; i < Particle_Number_PBCs; i++)
    {
      PBC = &(PBCs[i]);
      switch(PBC->type)
	{
	case PBC_TARGET:
	  fp = fopen(construct_filename(PBC->string_data), "w");
	  fclose(fp);
	  break;
	case PBC_IMPERMEABLE:
	  /* This is is sort of done "in reverse", so we just note it
	   * and catch it later. */
	  EH(-1, "CGM not supported, IMPERMEABLE PBC.");
	  /*
	  fprintf(stderr, "PBC %d is IMPERMEABLE, Num_Impermeable_PBCs = %d\n", i, Num_Impermeable_PBCs);
	  */
	  break;
	case PBC_OUTFLOW:
	case PBC_SOURCE:
	case PBC_FREESTREAM_SOURCE:
	  /* Nothing to do... */
	  break;
	default:
	  sprintf(err_msg, "Unknown PBC type (%d).", PBC->type);
	  EH(-1, err_msg);
	}
    }

}


/* The set of routines:
 *   load_element_information()
 *   load_element_node_coordinates()
 *   load_field_variables_at_xi()
 *
 * form a hierarchy of calls that depend on each other.  Through one
 * call, to the function that you "mean", you should get all of the
 * other structures properly filled in. */

/* See comment above.  The static global variables associated with
 * this routine are last_element_loaded and velo_interp. */
static void
load_element_information(const int elem_id)
{
  int err;
  int i;

  if(last_element_loaded == elem_id)
    return;
  err = load_elem_dofptr(elem_id, static_exo, static_x, static_x_old,
                         static_xdot, static_xdot_old, 0);
  EH(err, "Error from load_elem_dof_ptr()");
  
  /* Not sure if this call is necessary ... I don't know exactly what
   * element information I use laster on. */
  ei[pg->imtrx]->ielem = elem_id;
  load_ei(ei[pg->imtrx]->ielem, static_exo, 0, pg->imtrx);

  velo_interp = -1;
  for(i = 0; i < Num_Basis_Functions; i++)
    if(pd_glob[ei[pg->imtrx]->mn]->i[pg->imtrx][VELOCITY1] == bfd[i]->interpolation) velo_interp = i;
  if(velo_interp == -1)
    EH(-1, "Could not find velocity interpolation function.");

  last_element_loaded = elem_id;
}


/* See comment above.  This routine loads the corner point coordinates
 * for the requested element.  The static global variables associated
 * with this routine are last_elements_nodes_loaded and
 * node_coord[][]. */
static void
load_element_node_coordinates(const int elem_id)
{
  int i, j, k;
  int idof, node_index, var;

  if(elem_id == last_elements_nodes_loaded)
    return;
  load_element_information(elem_id);

  i = Proc_Connect_Ptr[elem_id];
  for(j = 0; j < mdim; j++)
    {
      if(pd->e[pg->imtrx][R_MESH1])
	{
	  var = MESH_DISPLACEMENT1 + j;
	  for(k = 0; k < nodes_per_element; k++)
	    {
	      idof = ei[pg->imtrx]->ln_to_dof[var][k];
	      node_index = Proc_Elem_Connect[i + k];
	      node_coord[k][j] = Coor[j][node_index] + *esp->d[j][idof];
	    }
	}
      else
	for(k = 0; k < nodes_per_element; k++)
	  {
	    node_index = Proc_Elem_Connect[i + k];
	    node_coord[k][j] = Coor[j][node_index];
	  }
    }
  last_elements_nodes_loaded = elem_id;
}


/* This routine filles an array of node indices for a particular side
 * of a hex element.  It does *not* guarantee a particular
 * orientation(counter-/clockwise), but it has one of them.  It also
 * guarantees that nodes 0 and 2 are opposite, as are 1 and 3, which
 * is all I needed for stuff here... */ 
static void
fill_hex_side_indices(const int side_id,
		      int indices[4])
{
  switch(side_id)
    {
    case 0: indices[0] = 0; indices[1] = 1; indices[2] = 5; indices[3] = 4;
      break;
    case 1: indices[0] = 1; indices[1] = 2; indices[2] = 6; indices[3] = 5;
      break;
    case 2: indices[0] = 2; indices[1] = 3; indices[2] = 7; indices[3] = 6;
      break;
    case 3: indices[0] = 0; indices[1] = 3; indices[2] = 7; indices[3] = 4;
      break;
    case 4: indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
      break;
    case 5: indices[0] = 4; indices[1] = 5; indices[2] = 6; indices[3] = 7;
      break;
    default:
      EH(-1, "Bad side_id.");
      break;
    }
}

/* See comment above.  This routine will pick a uniformly distributed
 * location from the requested element and side.
 *
 * NOTE: this routine assumes the side/face is linear.  It would be a
 * big pain to do this for actual curved elements, etc...
 */
static void
select_random_side_location(const int elem_id,
			    const int side_id,
			    dbl * x,
			    dbl * xi)
{
  int i, j, indices[4];
  dbl r, r2, x_edge[2][DIM];

  /* initialize indices */
  for (i = 0; i < 4; i++) {
    indices[i] = 0;
  }

  load_element_node_coordinates(elem_id);
  
  if(mdim == 2)
    {
      r = drand48();
      for(i = 0; i < mdim; i++)
	{
	  x[i] = node_coord[side_id][i];
	  x[i] += r * (node_coord[(side_id+1)%4][i] - x[i]);
	}
    }
  else
    {
      fill_hex_side_indices(side_id, indices);
      r = drand48();
      r2 = drand48();
      for(j = 0; j < DIM; j++)
	{
	  x_edge[0][j] = node_coord[indices[0]][j];
	  x_edge[1][j] = node_coord[indices[2]][j];
	  x_edge[0][j] += r * (node_coord[indices[1]][j] - x_edge[0][j]);
	  x_edge[1][j] += (1.0 - r) * (node_coord[indices[3]][j] - x_edge[1][j]);
	  x[j] = x_edge[0][j];
	  x[j] += r2 * (x_edge[1][j] - x[j]);
	}
      /*
      fprintf(stderr, "CREATED particle at (%g,%g,%g)\n",
	      x[0], x[1], x[2]);
      */
    }

  /* Fill in the element coordinates. */
  get_element_xi_newton(elem_id, x, xi);
}


/* This routine returns the area (3D) or length (2D) of an element's
 * side.
 *
 * NOTE: although it returns the side length for CYLINDRICAL, one
 * should modify any random selection routines to account for changes
 * in radii.  I haven't done this, yet. */
static dbl
get_side_area(const int elem_id,
	      const int side_id)
{
  int i, indices[4];
  dbl d1[DIM], d2[DIM], dcross[DIM], dnorm;

  /* initialize indices */
  for (i = 0; i < 4; i++) {
    indices[i] = 0;
  }

  load_element_node_coordinates(elem_id);
  if(mdim == 2)
    return my_distance(&(node_coord[side_id][0]),
		       &(node_coord[(side_id+1)%4][0]),
		       mdim);
  else
    {
      fill_hex_side_indices(side_id, indices);
      /* Get diagonals... */
      for(i = 0; i < 3; i++)
	{
	  d1[i] = node_coord[indices[2]][i] - node_coord[indices[0]][i];
	  d2[i] = node_coord[indices[3]][i] - node_coord[indices[1]][i];
	}
      cross_really_simple_vectors(d1, d2, dcross);
      dnorm = nnorm(DIM, dcross);

      /*
      printf("element %d, side %d, has area = %g\n",
	     elem_id, side_id, 0.5*dnorm);
      */
      
      return 0.5 * dnorm;
    }
}


/* See comment above.  This routine will load the field variables at a
 * passed element/xi location.  The static global variable(s)
 * associated with this routine are last_elements_fv_loaded and
 * last_fvs_xi_loaded.
 *
 * This routine, or one below it, will one day have reference to
 * Particle_Model to determine which field variables we really need to
 * load... */
static void
load_field_variables_at_xi(const int elem_id,
			   dbl * const xi)
{
  int err;
  int i, xis_equal;

  xis_equal = 1;
  for(i = 0; i < mdim; i++)
    if(xi[i] != last_fvs_xi_loaded[i])
      xis_equal = 0;
  if(last_elements_fv_loaded == elem_id && xis_equal)
    return;

  if(last_elements_fv_loaded != elem_id)
    load_element_information(elem_id);

  err = load_basis_functions(xi, bfd);
  EH(err, "Error from load_basis_functions()");

  err = beer_belly();
  EH(err, "Error from beer_belly()");

  err = load_fv();
  EH(err, "Error from load_fv()");
  if(TimeIntegration == STEADY)
    {
      /* The fv_old structures don't even get filled if we're a STEADY run. */
      if(pd->v[pg->imtrx][VELOCITY1])
	memcpy(fv_old->v, fv->v, DIM * sizeof(dbl));
      if(pd->v[pg->imtrx][ENORM])
	fv_old->Enorm = fv->Enorm;
    }

  /* When I put in deformable meshes, I need to include these two
   * calls, too. */

  err = load_bf_grad();
  EH(err, "Error from load_bf_grad()");

  err = load_fv_grads();
  EH(err, "Error from load_fv_grads()");

  last_elements_fv_loaded = elem_id;
  for(i = 0; i < mdim; i++)
    last_fvs_xi_loaded[i] = xi[i];
}


/* This routine finds the element coordinates (xi) for a specified
 * user coordinate (x).  It works for both 2D and 3D.  If the
 * coordinate does not lie on this element, then xi will reflect that
 * (some xi coordinate will be < -1 or > 1.  This is similar to a
 * routine that Randy wrote invert_isoparametric_map(), and I need to
 * ask him what the rotational tensor is.  If convergence cannot occur
 * then a value of -1 is returned (and a bunch of status data is
 * printed), otherwise 0 is returned.
 */
static int
get_element_xi_newton(const int elem_id,
		      const dbl * const x,
		      dbl * xi)
{
  dbl update[DIM], resid[DIM], resid_norm, xi_norm;
  dbl termination_residual;
  const dbl tol = 1.0e-12;
  int i, j, si, iter;
  struct Basis_Functions *map_bf;
  
  si = in_list(pd->IntegrationMap, 0, Num_Interpolations, Unique_Interpolations);
  map_bf = bfd[si];

  load_field_variables_at_xi(elem_id, xi);
  memset(resid, 0, DIM * sizeof(dbl));
  for(i = 0; i < mdim; i++)
    resid[i] = fv->x[i] - x[i];
  resid_norm = nnorm(mdim, resid);

  xi_norm = nnorm(mdim, xi);
  termination_residual = MAX(tol, tol * xi_norm*xi_norm);


  /*
  fprintf(stderr, "  GEXN Initial values: ");
#ifdef PARALLEL
  fprintf(stderr, " Proc %d ", ProcID);
#endif
  fprintf(stderr, "x = (%g,%g,%g), xi = (%g,%g,%g), fv->x = (%g,%g,%g), resid = %g (%g,%g,%g)\n",
	  x[0], x[1], x[2], xi[0], xi[1], xi[2], fv->x[0], fv->x[1], fv->x[2], resid_norm, resid[0], resid[1], resid[2]);
  fflush(stderr);
  */
  
  iter = 0;
  while(resid_norm > termination_residual &&
	iter < 100)
    {
      memset(update, 0, DIM * sizeof(dbl));
      for(i = 0; i < mdim; i++)
	for(j = 0; j < mdim; j++)
	  update[i] += map_bf->B[j][i] * resid[j];

      for(i = 0; i < mdim; i++)
	xi[i] -= update[i];

      load_field_variables_at_xi(elem_id, xi);
      for(i = 0; i < mdim; i++)
	resid[i] = fv->x[i] - x[i];
      resid_norm = nnorm(mdim, resid);

      /* If we're looking at a particle that's definitely out of this
       * element, account for the noise in locating it... */
      xi_norm = nnorm(mdim, xi);
      termination_residual = MAX(tol, tol * xi_norm*xi_norm);
      iter++;
      /*
      fprintf(stderr, "  GEXN update = (%g,%g,%g), xi = (%g,%g,%g), resid = (%g,%g,%g), resid norm = %g\n",
	     update[0], update[1], update[2], xi[0], xi[1], xi[2], resid[0], resid[1], resid[2], resid_norm);
      fprintf(stderr, "  GEXN fv->x = (%g,%g,%g)\n", fv->x[0], fv->x[1], fv->x[2]);
      fprintf(stderr, "  GEXN termination_residual = %g\n", termination_residual);
*/
    }

  if(iter == 100 && resid_norm > termination_residual)
    {
      fprintf(stderr, "\nget_element_xi_newton() says:\n");
      fprintf(stderr, "Searching from element %d for point (%g,%g,%g).\n", elem_id, x[0], x[1], x[2]);
      fprintf(stderr, "last update = (%g,%g,%g), xi = (%g,%g,%g)\nresid = (%g,%g,%g), resid norm = %g\n",
	     update[0], update[1], update[2], xi[0], xi[1], xi[2], resid[0], resid[1], resid[2], resid_norm);
      fprintf(stderr, "fv->x = (%g,%g,%g), fv->v = (%g,%g,%g)\n", fv->x[0], fv->x[1], fv->x[2], fv->v[0], fv->v[1], fv->v[2]);
      return -1;
    }

  /*
  fprintf(stderr, "  Leaving gexn() with x = (%g,%g,%g), xi = (%g,%g,%g)\n", fv->x[0], fv->x[1], fv->x[2], xi[0], xi[1], xi[2]);
  fprintf(stderr, "  GEXN terminated with resid = %g, steps = %d\n\n",
	 resid_norm, iter);
  */

  return 0;
}


/* This routine returns the minimum side length for the requested
 * element.  Currently used for particle CFL conditions.  There is an
 * obvious speed up for static meshes where we don't recompute these
 * values every time we're called, just pre-compute on a per-element
 * basis.
 */
static dbl
get_element_minimum_side_length(const int elem_id)
{
  int edge_id;
  dbl len = 1.0e+10;

  load_element_node_coordinates(elem_id);
  if(mdim == 2)
    {
      for(edge_id = 0; edge_id < 4; edge_id++)
	len = MIN(len, get_side_area(elem_id, edge_id));
    }
  else
    {
      len = MIN(len, my_distance(&(node_coord[0][0]), &(node_coord[1][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[0][0]), &(node_coord[3][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[0][0]), &(node_coord[4][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[1][0]), &(node_coord[2][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[1][0]), &(node_coord[5][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[2][0]), &(node_coord[3][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[2][0]), &(node_coord[6][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[3][0]), &(node_coord[7][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[4][0]), &(node_coord[5][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[4][0]), &(node_coord[7][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[5][0]), &(node_coord[6][0]), mdim));
      len = MIN(len, my_distance(&(node_coord[6][0]), &(node_coord[7][0]), mdim));
    }
  return len;
}
			    


/* This routine writes TECPLOT zone info.  This gets called at the
 * beginning of every compute_particle() step if we're using the
 * TECPLOT output format.  The intention here is that the user is
 * either (1) outputting particles at the same output times as the
 * regular Goma output (i.e., Output stride = 1), or (2) Still
 * essentially streaming the particles to file, but they are in
 * separate zones (and may have multiple particle time steps within a
 * single zone).  If we're computing after a steady-state calculation,
 * then the Output Time Step is used.
 *
 * The logic for when we compute this is confusing.  I'm worrying
 * primarly about the Output Stride = 1 case and the post-steady-state
 * Output Time Step case.
 */
static void
output_TECPLOT_zone_info(const dbl end_time,
			 const int goma_time_step,
			 const int force)
{
  int i, particles_exist;

  /* This isn't exactly right, but we cannot predict if we'll have
   * particles at the end of this step or not.  So assume we'll generally
   * always have them, or always *not* have them.
   */
#ifdef PARALLEL
  particles_exist = num_particles > 0;
  /*
  particles_exist = local_num_particles > 0;
  */
#else
  particles_exist = num_particles;
#endif

  if(!Particle_Number_Sample_Types)
    return;

  if(!particles_exist && !force)
    return;

  if((TimeIntegration == TRANSIENT && Particle_Output_Stride > 0 && (!((goma_time_step+1) % Particle_Output_Stride))) ||
     (TimeIntegration == STEADY && Particle_Output_Time_Step > 0.0))
    for(i = 0; i < Particle_Number_Sample_Types; i++)
      fprintf(pa_fp[i], "ZONE T=\"%g\"\n", end_time);
}


/* This routine handles quite a lot.  There are two basic kinds of
 * output to be performed: stride-based, and time-based.  For
 * stride-based, the Goma successful timestep needs to be known.  For
 * time-based the current particle time and last particle time need to
 * be known.  Furthermore, there are two kinds of outputs to be
 * considered.  The first is a sample type, which possibly has
 * associated primitive variables.  The second is a global type which
 * should just be a full location output for all particles.  Whew.
 *
 * You may be asking yourself why we do this much work for every
 * particle?  Well, for low-particle-high-work systems the particles
 * will have their own variable timestepping, etc., and this is the
 * better way to do it.  I probably should implement another whole
 * output method where the particle loop is the internal one for more
 * "DSMC"-ish simulations with many many particles.
 */
static void
output_a_particle(particle_t * const p,
		const dbl last_particle_time,
		const dbl particle_time,
		const int goma_time_step,
		const int last_step,
		const int force_output)
{
  int i, output_sample_number;
  int output_n, start_n, end_n;
  int strided_output, last_goma_step;
  dbl time_fraction, output_time;
  dbl x_time[DIM], real_data_time;

  
  /* The nice thing about strided output is that there is no
   * linearization by time fractions required... */
  if(Particle_Output_Stride > 0)
    {
      if((!((goma_time_step+1) % Particle_Output_Stride) ||
	  force_output)
	 && last_step)
	{
	  if((output_sample_number = p->output_sample_number) != -1)
	    {
	      for(i = 0; i < pdim; i++)
		fprintf(pa_fp[output_sample_number], " %12g ", p->x_old[i]);
	      fprintf(pa_fp[output_sample_number], " %12g", particle_time);
	      for(i = 0; i < Particle_Number_Output_Variables[output_sample_number]; i++)
		fprintf(pa_fp[output_sample_number], " %12g", p->real_data[i]);
#ifdef PARALLEL
	      fprintf(pa_fp[output_sample_number], " %d", p->owning_proc_id);
#endif
	      fprintf(pa_fp[output_sample_number], "\n");
	    }
	}
    }
  else if((output_sample_number = p->output_sample_number) != -1 ||
	  force_output)
    {
      start_n = (int)floor(last_particle_time/Particle_Output_Time_Step);
      end_n = (int)floor(particle_time/Particle_Output_Time_Step);
      if(TimeIntegration == STEADY) /* Need to take care of the "off by one" errors that can happen with the
				       floating point operations directly above, for post-steady calculations. */
	end_n = start_n = goma_time_step;
      if(end_n == start_n && force_output) end_n++;
      for(output_n = start_n; output_n < end_n; output_n++)
	{
	  if(particle_time == 0.0)
	    {
	      output_time = 0.0;
	      time_fraction = 1.0;
	    }
	  else
	    {
	      output_time = (output_n + 1) * Particle_Output_Time_Step;
	      time_fraction = (output_time - last_particle_time)/
		(particle_time - last_particle_time);
	    }
	  for(i = 0; i < pdim; i++)
	    x_time[i] = (1.0-time_fraction) * p->x_old[i] + time_fraction * p->x[i]; 
	  if(output_sample_number != -1)
	    {
	      for(i = 0; i < pdim; i++)
		fprintf(pa_fp[output_sample_number], " %12g", x_time[i]);
	      fprintf(pa_fp[output_sample_number], " %12g", output_time);
	      for(i = 0; i < Particle_Number_Output_Variables[output_sample_number]; i++)
		{
		  real_data_time = (1.0-time_fraction) * p->real_data_old[i] + time_fraction * p->real_data[i];
		  fprintf(pa_fp[output_sample_number], " %12g", real_data_time);
		}
#ifdef PARALLEL
	      fprintf(pa_fp[output_sample_number], " %d", p->owning_proc_id);
#endif
	      fprintf(pa_fp[output_sample_number], "\n");
	    }
	}
    }

  /* This is for the full output, which is always strided in terms of
   * either the post steady-state timesteps or the (transient) Goma
   * time step. */
  if(Particle_Full_Output_Stride > 0)
    {
      strided_output = !( (goma_time_step + 1) % Particle_Full_Output_Stride );
      last_goma_step = ((TimeIntegration == STEADY) && (goma_time_step  == Particle_Max_Time_Steps - 1)) ||
        ((TimeIntegration == TRANSIENT) && ((goma_time_step == tran->MaxTimeSteps - 1) || (particle_time == tran->TimeMax)));

      /*
      fprintf(stderr, "Proc%d: In oap/full, strided_output = %d, last_step = %d, last_goma_step = %d, force_output=%d\n", ProcID, strided_output, last_step, last_goma_step, force_output);
      */
      
      if((strided_output && last_step) ||
	 (last_goma_step && last_step) ||
	 force_output)
	{
	  for(i = 0; i < pdim; i++)
	    fprintf(pa_full_fp, " %12g", p->x[i]);
	  fprintf(pa_full_fp, " %12g", particle_time);
	  if(Particle_Model == SWIMMER_EXPLICIT || Particle_Model == SWIMMER_IMPLICIT)
	    {
	      fprintf(pa_full_fp, " %g", p->theta);
	      if(pdim == 3)
		fprintf(pa_full_fp, " %g", p->phi);
	    }
	  fprintf(pa_full_fp, " %d", p->owning_elem_id);
#ifdef PARALLEL
	  fprintf(pa_full_fp, " %d", p->owning_proc_id);
#endif
	  fprintf(pa_full_fp, "\n");
	}
    }
}


/* This routine stores the field variables (or other output
 * quantities) we're interested in on the particle.  It assumes that
 * load_fv, etc., as already been run.  I will change this when I
 * change the "load_fv" calls to just compute velocity. */
static void
store_particle_data(particle_t * const p)
{
  char *output_var;
  int i, j, output_sample_number;

  output_sample_number = p->output_sample_number;
  j = 0;
  for(i = 0; i < Particle_Number_Output_Variables[output_sample_number]; i++)
    {
      if(j == MAX_DATA_REAL_VALUES)
	dump1(EXIT, p, "Uh-oh, was about to overwrite the p->real_data array!");
      output_var = Particle_Output_Variables[output_sample_number][i];

      /* The reason I do this in such a speed-poor way
       * (because of all the strncmp's) is because I may want
       * to have totally user-defined output "values" at some
       * point.  And searching through the variable lists to
       * match the variables string name (i.e., "VX") would be
       * just as slow as this.  Ah well.*/

      /* Beware hiding longer-lengthed variables names. */
      /* Look into bc_colloc.c:load_variable() */
      if(!strncmp("T", output_var, 1))
	p->real_data[j++] = fv->T;
      else if(!strncmp("VX", output_var, 2))
	p->real_data[j++] = fv->v[0];
      else if(!strncmp("VY", output_var, 2))
	p->real_data[j++] = fv->v[1];
      else if(!strncmp("VZ", output_var, 2))
	p->real_data[j++] = fv->v[2];
      else
	dump2(EXIT, p, "Unknown particle output variable: %s.\n", output_var);
    }
}


/* This routine will compute a particle timestep size for this
 * particle.  It is bound above by the full Goma timestep.  The means
 * for determining this timestep size is entirely model dependent.
 * You better know what you're doing here... */
static dbl
compute_particle_dt(particle_t * const p,
		    const dbl max_dt)
{
  dbl dt = 0, tmp_dt = 0, particle_speed, fluid_speed, max_speed, min_side_length, mass, volume;
  /*  dbl fluid_mass; */
  dbl stokes_force_coeff, coulombic_force_coeff, gravity_force_coeff;
  dbl Re_p, Re_p_correction;
  dbl vel_diff, vel_rel[DIM];
  dbl forces[DIM], max_force;
  dbl coeff_a, coeff_b, coeff_c, coeff_d, CM_fact, dielectrophoretic_force_coeff;
  dbl force;
  int i, gravity_index;
  const dbl Particle_CFL_Maximum = 0.2;

  memset(forces, 0, sizeof(double)*DIM);

  /* Get the index into acceleration vectors for gravity. */
  if(pd->CoordinateSystem == CARTESIAN)
    gravity_index = 2*pdim - 1; /* y or z depending on 2D or 3D */
  else
    gravity_index = pdim; /* z in (z,r) */

  min_side_length = get_element_minimum_side_length(p->owning_elem_id);
  particle_speed = nnorm(DIM, p->v);
  fluid_speed = nnorm(DIM, fv->v);
  max_speed = MAX(particle_speed, fluid_speed);
  switch(Particle_Model)
    {
    case TRACER_EXPLICIT:
    case TRACER_IMPLICIT:
    case SWIMMER_EXPLICIT:
    case SWIMMER_IMPLICIT:
      if(max_speed == 0.0)
	dt = max_dt;
      else
	dt = Particle_CFL_Maximum * min_side_length / max_speed;
      dt = MIN(dt, max_dt);

      break;

    case INERTIAL_TRACER_EXPLICIT:
      EH(-1, "INERTIAL_TRACER_EXPLICIT offline -- need a better mover.");
      break;

    case INERTIAL_TRACER_IMPLICIT:
      /* These need to match the same values used in move_particle(). */
      volume = 4.0/3.0 * M_PIE * Particle_Radius * Particle_Radius * Particle_Radius;
      mass = Particle_Density * volume;
      /*      fluid_mass = mp->density * volume; */
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      gravity_force_coeff = Particle_Model_Data[1];

      for(i = 0; i < pdim; i++)
	vel_rel[i] = fv->v[i] - p->v[i];
      vel_diff = nnorm(pdim, vel_rel);
      Re_p = Particle_Density * vel_diff * 2.0 * Particle_Radius / mp->density;
      if(Re_p <= 0.1)
	Re_p_correction = 1.0 + 3.0/16.0 * Re_p;
      else
	Re_p_correction = 1.0 + 0.0565 * pow(Re_p, 0.525);

      if(Re_p > 500.0)
	fprintf(stderr, "Applying possibly improper Re_p_correction, Re_p = %g\n", Re_p);
      stokes_force_coeff *= Re_p_correction;

      /* Starting timestep based on particle info. */
      if(max_speed == 0.0)
	dt = 1.0e+10;
      else
	dt = Particle_CFL_Maximum * min_side_length / max_speed;

      /* Get the acceleration values and from that adjust timestep if necessary.  We assume that the fluid velocity
	 doesn't change appreciably over the timestep for CFL purposes. */
      for(i = 0; i < pdim; i++)
	forces[i] = stokes_force_coeff * (fv->v[i] - p->v[i])
	  - fv->grad_P[i] * volume;
      forces[gravity_index] += gravity_force_coeff * mass;

      max_force = fabs(forces[0]);
      for(i = 1; i < pdim; i++)
	if(fabs(forces[i]) > max_force)
	  max_force = fabs(forces[i]);

      if(max_force > 0.0)
	tmp_dt = sqrt(Particle_CFL_Maximum * min_side_length / max_force);
      else
	tmp_dt = dt;

      if(tmp_dt < dt / 10.0)
	fprintf(stderr, "Forces were dominant in choosing dt (old=%g, new=%g)\n", dt, tmp_dt);
      dt = MIN(dt, tmp_dt);

      /* Check for large forces due to fluid/particle velocity differences. */
/*       max_vel_diff =  0.0; */
/*       for(i = 0; i < pdim; i++) */
/* 	max_vel_diff = MAX(max_vel_diff, MAX(fabs(fv->v[i] - p->v[i]), fabs(fv_old->v[i] - p->v[i]))); */
/*       if(max_vel_diff > 0.0) */
/* 	tmp_dt = sqrt(Particle_CFL_Maximum * min_side_length / (stokes_force_coeff / mass * max_vel_diff)); */
/*       else */
/* 	tmp_dt = 1.0e+10; */
/*       dt = MIN(dt, tmp_dt); */

      /* Check for large gravity forces. */
/*       if(gravity_force_coeff != 0.0) */
/* 	{ */
/* 	  relative_mass = MAX(0.0, (mass - fluid_mass) / mass); */
/* 	  if(relative_mass != 0.0) */
/* 	    tmp_dt = sqrt(Particle_CFL_Maximum * min_side_length / (gravity_force_coeff * relative_mass)); */
/* 	  else */
/* 	    tmp_dt = 1.0e+10; */
/* 	} */
/*       else */
/* 	tmp_dt = 1.0e+10; */
/*       if(tmp_dt < dt / 10.0) */
/* 	fprintf(stderr, "Gravity was dominant in choosing dt (old=%g, new=%g)\n", dt, tmp_dt); */
/*       dt = MIN(dt, tmp_dt); */

      break;

    case CHARGED_TRACER_EXPLICIT:
      EH(-1, "CHARGED_TRACER_EXPLICIT offline -- need a better mover.");
      break;

    case CHARGED_TRACER_IMPLICIT:
      /* These need to match the same values used in move_particle(). */
      mass = 4.0/3.0 * M_PIE * Particle_Density * Particle_Radius * Particle_Radius * Particle_Radius;
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      coulombic_force_coeff = -Particle_Model_Data[1];

      /* Starting timestep based on particle info. */
      if(max_speed == 0.0)
	dt = 1.0e+10;
      else
	dt = Particle_CFL_Maximum * min_side_length / max_speed;

      /* Get the current grad_V value. */
      /*
      fill_my_fv_old();
      */
      
      /* Check for large forces due to fluid/particle velocity differences. */
      /* This needs to be modified to use old values of v and V, too. */
      max_force = 0.0;
      for(i = 0; i < pdim; i++)
	{
	  force = fabs(stokes_force_coeff/mass * (fv->v[i] - p->v[i]) +
		       coulombic_force_coeff/mass * fv->grad_V[i]);
	  max_force = MAX(max_force, force);
	}
      if(max_force > 0.0)
	tmp_dt = sqrt(Particle_CFL_Maximum * min_side_length / max_force);
      dt = MIN(dt, tmp_dt);

      break;

    case DIELECTROPHORETIC_TRACER_IMPLICIT:

      /* The 6 coefficients for DIELECTROPHORETIC_TRACER_IMPLICIT are:
       * 1) Brownian motion (currently neglected, and should be zero)
       * 2) Electrical permittivity of the particle.
       * 3) Electrical permittivity of the medium (fluid).
       * 4) Electrical conductivity of the particle.
       * 5) Electrical conductivity of the medium (fluid).
       * 6) Frequency of the AC field.
       * 7) Volt conversion factor.
       */

      /* These need to match the same values used in move_particle(). */
      mass = 4.0/3.0 * M_PIE * Particle_Density * Particle_Radius * Particle_Radius * Particle_Radius;
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      /* Compute the Clausius-Mossoti factor, Re(K(omega)) */
      coeff_a = Particle_Model_Data[1] - Particle_Model_Data[2];
      coeff_b = (Particle_Model_Data[4] - Particle_Model_Data[3]) / Particle_Model_Data[5];
      coeff_c = Particle_Model_Data[1] + 2.0 * Particle_Model_Data[2];
      coeff_d = (Particle_Model_Data[3] + 2.0 * Particle_Model_Data[4]) / Particle_Model_Data[5];
      CM_fact = (coeff_a * coeff_c - coeff_b * coeff_d) / (coeff_c * coeff_c + coeff_d * coeff_d);
      dielectrophoretic_force_coeff = 2.0 * M_PIE * Particle_Radius * Particle_Radius * Particle_Radius	*
	Particle_Model_Data[2] * CM_fact * Particle_Model_Data[6] * Particle_Model_Data[6];

      /*
      fprintf(stderr, "New dt calculation\n");
      fprintf(stderr, "\tCM_fact = %g, dfc=%g\n", CM_fact, dielectrophoretic_force_coeff);
      */
      
      /* Starting timestep based on particle info. */
      if(max_speed == 0.0)
	dt = 1.0e+10;
      else
	dt = Particle_CFL_Maximum * min_side_length / max_speed;

      /*
      fprintf(stderr, "\tloc = (%g,%g,%g)\n", fv->x[0], fv->x[1], fv->x[2]);
      fprintf(stderr, "\tCFL dt = %g\n", dt);
      */      

      /* Get the current grad_V value. */
      /*
      fill_my_fv_old();
      */
      
      /* Check for large forces due to fluid/particle velocity differences. */
      /* This needs to be modified to use old values of v and V, too. */
      max_force = 0.0;
      for(i = 0; i < pdim; i++)
	{
	  force = fabs(stokes_force_coeff/mass * (fv->v[i] - p->v[i]) +
		       dielectrophoretic_force_coeff/mass * 2.0 * fv->Enorm * fv->grad_Enorm[i]);

	  /*
	  fprintf(stderr, "\ti=%d sfc/m=%g sfc/m(%g)=%g\n", i, stokes_force_coeff/mass, fv->v[i] - p->v[i], stokes_force_coeff/mass*(fv->v[i]-p->v[i]));
	  fprintf(stderr, "\ti=%d dfc/m=%g dfc/m(2)(%g) = dfc/m(2)(%g)(%g)=%g\n", i, dielectrophoretic_force_coeff/mass, fv->Enorm * fv->grad_Enorm[i], fv->Enorm, fv->grad_Enorm[i], dielectrophoretic_force_coeff/mass * 2.0 * fv->Enorm * fv->grad_Enorm[i]);
	  */
	  
	  max_force = MAX(max_force, force);

	  /*
	  fprintf(stderr, "\ti=%d: f=%g, sf=%g, df=%g\n", i, force, stokes_force_coeff/mass*(fv->v[i]-p->v[i]), dielectrophoretic_force_coeff/mass * 2.0 * fv->Enorm * fv->grad_Enorm[i]);
	  */
	}
      if(max_force > 0.0)
	{
	  tmp_dt = sqrt(Particle_CFL_Maximum * min_side_length / max_force);
	  /*
	  fprintf(stderr, "\tDielectrophoretic dt = %g\n", tmp_dt);
	  */
	}
      else
	tmp_dt = dt;
      dt = MIN(dt, tmp_dt);
      /*
      fprintf(stderr, "\tFinal dt = %g\n", dt);
      */
      break;

    default:
      dump1(EXIT, p, "Unknown Particle_Model.");
    }

  dt = MIN(dt, max_dt);
  return dt;
}


/* This routine computes the particle trajectory for the current
 * particle timestep.  There are explicit moves, implicit moves,
 * different forces to consider, etc.  The routine returns the number
 * of newton iterations required. */
static int
move_particle(particle_t * const p,
	      dbl * max_particle_dt,
	      dbl global_start_time,
	      dbl global_time_step)
{
  int i, j, return_value = -1;
  int gravity_index;
  dbl ext_velocity[DIM];	/* Additional velocity applied to particle. */
  dbl total_velocity[DIM];	/* total velocity at particle. */
  dbl time_fraction;
  dbl current_v[DIM];		/* For implicit methods, this is current velocity. */
  dbl vel_rel[DIM];		/* Relative velocity. */
  dbl current_grad_V[DIM];	/* For voltage field-based forces. */
  dbl current_Enorm;            /* Time-accurate |E| for dielectrophoresis. */
  dbl current_grad_Enorm[DIM];	/* Time-accurate grad(|E|) for dielectrophoresis. */
  int newton_iterations, newton_iterations_2;
  dbl J[6][6], b[6], x_sol[6], b_norm, relax_factor, update_norm;
  dbl rcubed;
  dbl volume, mass;
  /* dbl fluid_mass; */
  dbl weiner_angle, weiner_speed, beta;
  dbl Re_p, Re_p_correction, vel_diff;
  dbl stokes_force_coeff, coulombic_force_coeff, gravity_force_coeff;
  dbl coeff_a, coeff_b, coeff_c, coeff_d, CM_fact, dielectrophoretic_force_coeff;
  dbl J_2x2[2][2];
  
  dbl particle_dt, min_particle_dt;
  dbl xi_max_n, xi_max_np1, xi_target, dt_update = 0.0, d_xi_d_dt, xi_tolerance;;

  particle_dt = *max_particle_dt;
  min_particle_dt = MIN(1.0e-2, particle_dt);
  time_fraction = particle_dt / global_time_step;
  xi_tolerance = 1.0e-3;
  xi_target = 1.0 + xi_tolerance;

  load_field_variables_at_xi(p->owning_elem_id, p->xi);

  /* Everything needs the fluid velocity, so I can do this here.  This
   * only works for nondeformable meshes... */
  for(i = 0; i < DIM; i++)
    current_v[i] = (1.0 - time_fraction) * fv_old->v[i] + time_fraction * fv->v[i];

  rcubed = Particle_Radius * Particle_Radius * Particle_Radius;
  volume = 4.0/3.0 * M_PIE * rcubed;
  mass = volume * Particle_Density;

  /* Get the index into Newton RHS's for gravity. */
  if(pd->CoordinateSystem == CARTESIAN)
    gravity_index = 2*pdim - 1; /* y or z depending on 2D or 3D */
  else
    gravity_index = pdim; /* z in (z,r) */

  switch(Particle_Model)
    {
    case TRACER_EXPLICIT:
    case SWIMMER_EXPLICIT:
      memset(ext_velocity, 0, DIM * sizeof(dbl));
      /* Currently only SWIMMER_EXPLICIT has a real stochastic component... */
      if(Particle_Model == SWIMMER_EXPLICIT)
	{
	  /* 0th parameter is sigma for orientation angle in radians/second.
	   * 1st parameter is sigma (in fractions of swimming speed) for swimming speed in L/second.
	   * 2nd parameter is the swimming speed.
	   */
	  beta = 8.0 * M_PIE * mp->viscosity * rcubed;
	  /* Including moment of inertia => */
	  beta *= 2.0 / 5.0 * mass * Particle_Radius * Particle_Radius;
	  
	  weiner_angle = drand_truncated_normal(0.0, Particle_Model_Data[0] * sqrt(particle_dt));
	  weiner_speed = drand_truncated_normal(Particle_Model_Data[2],
						Particle_Model_Data[1] * Particle_Model_Data[2] * sqrt(particle_dt));

	  /*
	  fprintf(stderr, "Proc%d: beta = %g, mass = %g\n", ProcID, beta, mass);
	  fprintf(stderr, "Proc%d: torque = %g\n", ProcID, 980.0 * Particle_Radius * 0.05 * mass * sin(p->theta_old) / beta);
	  */
	  
	  if(mdim == 2)
	    {
	      /* For 2D, only one angle, theta.  theta = pi is up. */
	      p->theta_old = p->theta;
	      p->theta += particle_dt * (980.0 * Particle_Radius * 0.05 * mass * sin(p->theta_old) / beta) /* gravity-induced torque */
		+ weiner_angle; /* stochastic component */
	      /* Drop that ODE crap, just set it ... */
	      p->theta = M_PIE + weiner_angle;
	      p->theta = fmod(p->theta, 2.0 * M_PIE);
	      if(p->theta < 0.0) p->theta += 2.0 * M_PIE;
	      
	      ext_velocity[0] = weiner_speed * cos(p->theta - 0.5 * M_PIE);
	      ext_velocity[1] = weiner_speed * sin(p->theta - 0.5 * M_PIE);
	    }
	  else
	    {
	      p->theta_old = p->theta;
	      p->phi_old = p->phi;
	      
	      p->theta = drand48() * 2.0 * M_PIE;
	      p->phi = weiner_angle;

	      /* While I'm just setting it to (0, 0, 1) + stochastic,
	       * I don't care about the actual orientation for
	       * theta/phi, but I need to fix this when I actually go
	       * to a ODE form... */
	      ext_velocity[0] = weiner_speed * sin(p->phi) * cos(p->theta);
	      ext_velocity[1] = weiner_speed * sin(p->phi) * sin(p->theta);
	      ext_velocity[2] = weiner_speed * cos(p->phi);
	    }
	}
      else
	if(Particle_Model_Data[0] != 0.0)
	  EH(-1, "Set up random components!");

      memcpy(total_velocity, current_v, DIM * sizeof(dbl));
      for(i = 0; i < DIM; i++)
	total_velocity[i] += ext_velocity[i];

      if(Num_Impermeable_PBCs)
	enforce_impermeable_PBCs(p, total_velocity, particle_dt);

      for(i = 0; i < pdim; i++)
	{
	  p->x[i] += total_velocity[i] * particle_dt;
	  p->v[i] = total_velocity[i];
	}

      return_value = 0;
      break;

    case CHARGED_TRACER_EXPLICIT: /* This is *always* INERTIAL!! */
      /* Get the current grad_V value. */
      fill_my_fv_old();
      for(i = 0; i < pdim; i++)
	current_grad_V[i] = (1.0 - time_fraction) * (my_fv_old->grad_V[i]) + time_fraction * (fv->grad_V[i]);

      /* Fill the coefficients. */
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      coulombic_force_coeff = - Particle_Model_Data[1];
      
      for(i = 0; i < pd->Num_Dim; i++)
	total_velocity[i] = p->v[i] + particle_dt *
	  (stokes_force_coeff/mass * (current_v[i] - p->v[i]) +
	   coulombic_force_coeff/mass * current_grad_V[i]);

      if(Num_Impermeable_PBCs)
	enforce_impermeable_PBCs(p, total_velocity, particle_dt);

      for(i = 0; i < pd->Num_Dim; i++)
	{
	  p->x[i] += total_velocity[i] * particle_dt;
	  p->v[i] = total_velocity[i];
	}

      return_value = 0;
      break;

    case CHARGED_TRACER_IMPLICIT:
      /* Get the current grad_V value. */
      fill_my_fv_old();
      for(i = 0; i < pdim; i++)
	current_grad_V[i] = (1.0 - time_fraction) * (my_fv_old->grad_V[i]) + time_fraction * (fv->grad_V[i]);

      /* These need to match the same values used in compute_particle_dt(). */
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      coulombic_force_coeff = - Particle_Model_Data[1];

      for(i = 0; i < pdim; i++)
	{
	  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	  b[pdim+i] = p->v[i] - p->v_old[i] - particle_dt *
	    (stokes_force_coeff/mass * (current_v[i] - p->v[i]) +
	     coulombic_force_coeff/mass * current_grad_V[i]);
	}
      b_norm = nnorm(2*pdim, b);

      /*
      b[0] = p->x[0] - p->x_old[0] - particle_dt * p->v[0];
      b[1] = p->x[1] - p->x_old[1] - particle_dt * p->v[1];
      b[2] = p->v[0] - p->v_old[0] - particle_dt *
	(stokes_force_coeff/mass * (current_v[0] - p->v[0]) +
	 coulombic_force_coeff/mass * current_grad_V[0]);
      b[3] = p->v[1] - p->v_old[1] - particle_dt *
	(stokes_force_coeff/mass * (current_v[1] - p->v[1]) +
	 coulombic_force_coeff/mass * current_grad_V[1]);
      b_norm = nnorm(4, b);
      */
      
      update_norm = 1.0e+10;
      newton_iterations = 0;
      while(b_norm > 1.0e-8 &&
	    update_norm > 1.0e-8 &&
	    newton_iterations < 100)
	{
	  memset(J, 0, 36 * sizeof(dbl));
	      
	  /* ACK!  I need d(grad_V[i])/d(x[j]) !!  V = voltage
	   * Soooooooo, since the fv->grad_v only shows up in the
	   * Jacobian, do I really need to get current_grad_v (and
	   * thus create a my_fv_old->grad_v?), or will this still
	   * converge just fine?  I'll assume it converges just fine,
	   * but if there are problems, then replace the fv->grad_v
	   * values with current_grad_v values...
	   */
	  for(i = 0; i < pdim; i++)
	    {
	      J[i][i] = 1.0;
	      J[i][pdim+i] = -particle_dt;
	      J[pdim+i][pdim+i] = 1.0 + particle_dt * stokes_force_coeff/mass;
	      for(j = 0; j < pdim; j++)
		J[pdim+i][j] = -particle_dt * stokes_force_coeff/mass * fv->grad_v[j][i];
	    }

	  /*
	  fprintf(stderr, "J | b = % 12.6g % 12.6g % 12.6g % 12.6g | % 12.6g\n",
		  J[0][0], J[0][1], J[0][2], J[0][3], b[0]);
	  for(i = 1; i < 4; i++)
	    fprintf(stderr, "        % 12.6g % 12.6g % 12.6g % 12.6g | % 12.6g\n",
		    J[i][0], J[i][1], J[i][2], J[i][3], b[i]);
	  */
	  
	  solve_NxN_system(&J[0][0], b, x_sol, 2*pdim, 6);

	  /*
	  fprintf(stderr, "update = (%g,%g,%g,%g)\n",
		  -x_sol[0], -x_sol[1], -x_sol[2], -x_sol[3]);
	  */
	  
	  update_norm = nnorm(2*pdim, x_sol);
	  /* include min_side_length, too? */
	  relax_factor = 1.0 - (update_norm / (update_norm + nnorm(pdim, p->v)));
	  if(relax_factor < 0.1 && relax_factor != 0.0)
	    fprintf(stderr, "NEWTON Relax factor = %g\n", relax_factor);
	  if(relax_factor == 0.0)
	    relax_factor = 1.0;

	  for(i = 0; i < 2*pdim; i++)
	    x_sol[i] *= relax_factor;
	  for(i = 0; i < pdim; i++)
	    {
	      p->x[i] -= x_sol[i];
	      p->v[i] -= x_sol[pdim+i];
	    }

	  if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
	    dump1(EXIT, p, "Could not find element coordinates.");
	  /*
	  if(fabs(p->xi[0]) > 1.0 || fabs(p->xi[1]) > 1.0)
	    {
	      find_exit_wound(p->owning_elem_id, p, p->x_old, p->x, p->xi);
	    }
	  */

	  load_field_variables_at_xi(p->owning_elem_id, p->xi);
	  for(i = 0; i < pdim; i++)
	    current_v[i] = (1.0 - time_fraction) * (fv_old->v[i]) + time_fraction * (fv->v[i]);
	  fill_my_fv_old();
	  for(i = 0; i < pdim; i++)
	    current_grad_V[i] = (1.0 - time_fraction) * (my_fv_old->grad_V[i]) + time_fraction * (fv->grad_V[i]);

	  for(i = 0; i < pdim; i++)
	    {
	      b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	      b[pdim+i] = p->v[i] - p->v_old[i] - particle_dt *
		(stokes_force_coeff/mass * (current_v[i] - p->v[i]) +
		 coulombic_force_coeff/mass * current_grad_V[i]);
	    }
	  b_norm = nnorm(2*pdim, b);

	  newton_iterations++;
	}

      if(newton_iterations == 100 &&
	 b_norm > 1.0e-8 &&
	 update_norm > 1.0e-8)
	{
	  fprintf(stderr, "\nmove_particle() says:\n");
	  fprintf(stderr, "fv->x = (%g,%g,%g)\n",
		  fv->x[0], fv->x[1], fv->x[2]);
	  fprintf(stderr, "update = (%g,%g,%g,%g,%g,%g)\n", -x_sol[0], -x_sol[1], -x_sol[2], -x_sol[3], -x_sol[4], -x_sol[5]);
	  fprintf(stderr, "resid = (%g,%g,%g,%g,%g,%g), resid_norm = %g\n",
		  b[0], b[1], b[2], b[3], b[4], b[5], b_norm);
	  fprintf(stderr, "particle_dt = %g\n", particle_dt);
	  dump1(EXIT, p, "Failed to converge in 100 Newton iterations.");
	}

      return_value = newton_iterations;
      break;

    case DIELECTROPHORETIC_TRACER_IMPLICIT:
      /* Get the current grad_V value. */
      fill_my_fv_old();
      current_Enorm = (1.0 - time_fraction) * (fv_old->Enorm) + time_fraction * (fv->Enorm);
      for(i = 0; i < pdim; i++)
	current_grad_Enorm[i] = (1.0 - time_fraction) * (my_fv_old->grad_Enorm[i]) + time_fraction * (fv->grad_Enorm[i]);

      /* These need to match the same values used in compute_particle_dt(). */
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;

      /* Compute the Clausius-Mossoti factor, Re(K(omega)) */
      coeff_a = Particle_Model_Data[1] - Particle_Model_Data[2];
      coeff_b = (Particle_Model_Data[4] - Particle_Model_Data[3]) / Particle_Model_Data[5];
      coeff_c = Particle_Model_Data[1] + 2.0 * Particle_Model_Data[2];
      coeff_d = (Particle_Model_Data[3] + 2.0 * Particle_Model_Data[4]) / Particle_Model_Data[5];
      CM_fact = (coeff_a * coeff_c - coeff_b * coeff_d) / (coeff_c * coeff_c + coeff_d * coeff_d);
      dielectrophoretic_force_coeff = 2.0 * M_PIE * Particle_Radius * Particle_Radius * Particle_Radius	*
	Particle_Model_Data[2] * CM_fact * Particle_Model_Data[6] * Particle_Model_Data[6];

      for(i = 0; i < pdim; i++)
	{
	  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	  b[pdim+i] = mass * (p->v[i] - p->v_old[i]) - particle_dt *
	    (stokes_force_coeff * (current_v[i] - p->v[i]) +
	     dielectrophoretic_force_coeff * 2.0 * current_Enorm * current_grad_Enorm[i]);

	  /*
	    fprintf(stderr, "b[%d] = %g * (%g) - %g * ( %g * (%g) + %g * %g * %g)\n",
		  i, mass, p->v[i]-p->v_old[i], particle_dt, stokes_force_coeff, current_v[i] - p->v[i],
		  dielectrophoretic_force_coeff*2.0, current_Enorm, current_grad_Enorm[i]);
	  */
	}

      /* Force at least 1 iteration... */
      b_norm = 1.0e+10;
      xi_max_n = MAX(fabs(p->xi[0]), MAX(fabs(p->xi[1]), fabs(p->xi[2])));
      /*
      fprintf(stderr, "Initial xi_max = %.16g\n", xi_max_n);
      */
      
      /* vel is O(1e-3), dx is O(1e-6) */
      update_norm = 1.0e+10;
      newton_iterations = 0;
      while(b_norm > 1.0e-8 &&
	    update_norm > 1.0e-8 &&
	    newton_iterations < 100)
	{
	  /*
	  fprintf(stderr, "Beginning of newton iteration %d,\n", newton_iterations);
	  fprintf(stderr, "\tb_norm = %g (vs. 1.0e-8)\n", b_norm);
	  fprintf(stderr, "\tupdate_norm = %g (vs. 1.0e-8)\n", update_norm);
	  */
	  memset(J, 0, 36 * sizeof(dbl));

	  /* ACK!  I need d(grad_V[i])/d(x[j]) !!  V = voltage
	   * Soooooooo, since the fv->grad_v only shows up in the
	   * Jacobian, do I really need to get current_grad_v (and
	   * thus create a my_fv_old->grad_v?), or will this still
	   * converge just fine?  I'll assume it converges just fine,
	   * but if there are problems, then replace the fv->grad_v
	   * values with current_grad_v values...
	   */
	  for(i = 0; i < pdim; i++)
	    {
	      J[i][i] = 1.0;
	      J[i][pdim+i] = -particle_dt;
	      J[pdim+i][pdim+i] = mass + particle_dt * stokes_force_coeff;
	      for(j = 0; j < pdim; j++)
		J[pdim+i][j] = -particle_dt *
		  (stokes_force_coeff * fv->grad_v[j][i] +
		   dielectrophoretic_force_coeff * 2.0 * current_grad_Enorm[j] * current_grad_Enorm[i]);
	    }
	  /*
	  fprintf(stderr, "Going into solver, J=\n");
	  for(i = 0; i < 6; i++)
	    fprintf(stderr, "%12g %12g %12g %12g %12g %12g\n",
		    J[i][0], J[i][1], J[i][2], J[i][3], J[i][4], J[i][5]);
	  fprintf(stderr, "                   b=\n");
	  for(i = 0; i < 6; i++)
	    fprintf(stderr, "%12g\n", b[i]);
	  */
	  
	  solve_NxN_system(&J[0][0], b, x_sol, 2*pdim, 6);

	  update_norm = nnorm(2*pdim, x_sol);
	  /* include min_side_length, too? */
	  /* This doesn't work quite right when p->v -> 0, perhaps faster than update_norm... */

	  /*
	  relax_factor = 1.0 - (update_norm / (update_norm + nnorm(pdim, p->v)));
	  if(1)
	    fprintf(stderr, "NEWTON Relax factor = %g\n", relax_factor);
	  if(relax_factor == 0.0)
	    relax_factor = 1.0e-10;
	  */
	  relax_factor = 1.0;

	  /*
	  fprintf(stderr, "update before relaxation = (%g,%g,%g,%g,%g,%g)\n",
		  x_sol[0], x_sol[1], x_sol[2], x_sol[3], x_sol[4], x_sol[5]);
	  */
	  
	  for(i = 0; i < 2*pdim; i++)
	    x_sol[i] *= relax_factor;

	  /*
	  fprintf(stderr, "update after relaxation = (%g,%g,%g,%g,%g,%g)\n",
		  x_sol[0], x_sol[1], x_sol[2], x_sol[3], x_sol[4], x_sol[5]);
	  */
	  
	  for(i = 0; i < pdim; i++)
	    {
	      p->x[i] -= x_sol[i];
	      p->v[i] -= x_sol[pdim+i];
	    }

	  /*
	  fprintf(stderr, "new solution = (%g,%g,%g,%g,%g,%g)\n",
		  p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2]);
	  */
	  
	  if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
	    dump1(EXIT, p, "Could not find element coordinates.");

	  /*
	  fprintf(stderr, "new solution xi = (%g,%g,%g)\n",
		  p->xi[0], p->xi[1], p->xi[2]);
	  */
	  
	  load_field_variables_at_xi(p->owning_elem_id, p->xi);
	  for(i = 0; i < pdim; i++)
	    current_v[i] = (1.0 - time_fraction) * (fv_old->v[i]) + time_fraction * (fv->v[i]);
	  current_Enorm = (1.0 - time_fraction) * (fv_old->Enorm) + time_fraction * (fv->Enorm);
	  fill_my_fv_old();
	  for(i = 0; i < pdim; i++)
	    current_grad_Enorm[i] = (1.0 - time_fraction) * (my_fv_old->grad_Enorm[i]) + time_fraction * (fv->grad_Enorm[i]);

	  for(i = 0; i < pdim; i++)
	    {
	      b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	      b[pdim+i] = mass * (p->v[i] - p->v_old[i]) - particle_dt *
		(stokes_force_coeff * (current_v[i] - p->v[i]) +
		 dielectrophoretic_force_coeff * 2.0 * current_Enorm * current_grad_Enorm[i]);
	      /*
	      fprintf(stderr, "b[%d] = %g * (%g) - %g * ( %g * (%g) + %g * %g * %g\n",
		      i, mass, p->v[i]-p->v_old[i], particle_dt, stokes_force_coeff, current_v[i] - p->v[i],
		      dielectrophoretic_force_coeff*2.0, current_Enorm, current_grad_Enorm[i]);
	      */
	    }
	  b_norm = nnorm(2*pdim, b);

	  newton_iterations++;
	}
      /*
      fprintf(stderr, "Out of newton iteration %d,\n", newton_iterations);
      fprintf(stderr, "\tb_norm = %g (vs. 1.0e-8)\n", b_norm);
      fprintf(stderr, "\tupdate_norm = %g (vs. 1.0e-8)\n", update_norm);
      */
      
      xi_max_np1 = MAX(fabs(p->xi[0]), MAX(fabs(p->xi[1]), fabs(p->xi[2])));

      /* If we went out of the element, then solve the problem for the
       * proper dt so xi_max = 1.0 */
      newton_iterations_2 = 0;
      if(xi_max_np1 > xi_target)
	{
	  /*
	  fprintf(stderr, "XI_MAX xi_max_np1 = %.16g = ||(%12g,%12g,%12g)||_1\n", 
		  xi_max_np1, p->xi[0], p->xi[1], p->xi[2]);
	  */
	  b_norm = 1.0e+10;
	  update_norm = 1.0e+10;
	  while(((fabs(xi_max_np1 - xi_target) > xi_tolerance && /* not near the target */
		  particle_dt != min_particle_dt && /* used min step size twice in a row */
		  fabs(xi_max_np1 - xi_max_n) > xi_tolerance) || /* dt made no difference */
		 (b_norm >= 1.0e-8 &&
		  update_norm > 1.0e-8)) &&
		newton_iterations_2 < 100)
	    {
	      /*
	      fprintf(stderr, "XI_MAX xi_max_n = %.16g, xi_max_np1 = %.16g\n",
		      xi_max_n, xi_max_np1);
	      fprintf(stderr, "XI_MAX Changing PARTICLE_DT from %g to ", particle_dt);
	      */
	      if(newton_iterations_2 == 0)
		{
		  dt_update = particle_dt * (xi_target - xi_max_np1) / (xi_max_np1 - xi_max_n);
		}
	      else
		{
		  d_xi_d_dt = (xi_max_np1 - xi_max_n) / dt_update;
		  dt_update = (xi_target - xi_max_np1) / d_xi_d_dt;
		}
	      particle_dt += dt_update;
	      particle_dt = MAX(particle_dt, min_particle_dt);
	      time_fraction = particle_dt / global_time_step;

	      /*
	      fprintf(stderr, "%g (delta=%g), |xi_max_np1-xi_target| was %.16g\n",
		      particle_dt, dt_update, fabs(xi_target-xi_max_np1));
	      if(newton_iterations > 0)
		fprintf(stderr, "XI_MAX dt_update = %.16g, d_xi_d_dt = %.16g\n",
			dt_update, d_xi_d_dt);
	      */
	      
	      xi_max_n = xi_max_np1;

	      for(i = 0; i < pdim; i++)
		{
		  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
		  b[pdim+i] = mass * (p->v[i] - p->v_old[i]) - particle_dt *
		    (stokes_force_coeff * (current_v[i] - p->v[i]) +
		     dielectrophoretic_force_coeff * 2.0 * current_Enorm * current_grad_Enorm[i]);
		  /*		  
		  fprintf(stderr, "XI_MAX b[%d] = %g * (%g) - %g * ( %g * (%g) + %g * %g * %g)\n",
			  i, mass, p->v[i]-p->v_old[i], particle_dt, stokes_force_coeff, current_v[i] - p->v[i],
			  dielectrophoretic_force_coeff*2.0, current_Enorm, current_grad_Enorm[i]);
		  */
		}
	      memset(J, 0, 36 * sizeof(dbl));
	      for(i = 0; i < pdim; i++)
		{
		  J[i][i] = 1.0;
		  J[i][pdim+i] = -particle_dt;
		  J[pdim+i][pdim+i] = mass + particle_dt * stokes_force_coeff;
		  for(j = 0; j < pdim; j++)
		    J[pdim+i][j] = -particle_dt *
		      (stokes_force_coeff * fv->grad_v[j][i] +
		       dielectrophoretic_force_coeff * 2.0 * current_grad_Enorm[j] * current_grad_Enorm[i]);
		}
	      /*
	      fprintf(stderr, "XI_MAX Going into solver, J=\n");
	      for(i = 0; i < 6; i++)
		fprintf(stderr, "XI_MAX %12g %12g %12g %12g %12g %12g\n",
			J[i][0], J[i][1], J[i][2], J[i][3], J[i][4], J[i][5]);
	      fprintf(stderr, "XI_MAX                    b=\n");
	      for(i = 0; i < 6; i++)
		fprintf(stderr, "XI_MAX %12g\n", b[i]);
	      */	  
	      solve_NxN_system(&J[0][0], b, x_sol, 2*pdim, 6);

	      update_norm = nnorm(2*pdim, x_sol);

	      /* include min_side_length, too? */
	      /* This doesn't work quite right when p->v -> 0, perhaps faster than update_norm... */
	      /*
	      relax_factor = 1.0 - (update_norm / (update_norm + norm(pdim, p->v)));
	      if(1)
		fprintf(stderr, "XI_MAX NEWTON Relax factor = %g\n", relax_factor);
	      if(relax_factor == 0.0)
		relax_factor = 1.0e-10;
	      */
	      relax_factor = 1.0;

	      /*
		fprintf(stderr, "update before relaxation = (%g,%g,%g,%g,%g,%g)\n",
		x_sol[0], x_sol[1], x_sol[2], x_sol[3], x_sol[4], x_sol[5]);
	      */
	  
	      for(i = 0; i < 2*pdim; i++)
		x_sol[i] *= relax_factor;

	      /*
	      fprintf(stderr, "XI_MAX update after relaxation = (%g,%g,%g,%g,%g,%g)\n",
	      		      x_sol[0], x_sol[1], x_sol[2], x_sol[3], x_sol[4], x_sol[5]);
	      */
	      
	      for(i = 0; i < pdim; i++)
		{
		  p->x[i] -= x_sol[i];
		  p->v[i] -= x_sol[pdim+i];
		}

	      /*
	      fprintf(stderr, "XI_MAX new solution = (%g,%g,%g,%g,%g,%g)\n",
		      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2]);
	      */
	      if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
		dump1(EXIT, p, "XI_MAX Could not find element coordinates.");

	      xi_max_np1 = MAX(fabs(p->xi[0]), MAX(fabs(p->xi[1]), fabs(p->xi[2])));
	      /*
	      fprintf(stderr, "XI_MAX new solution |xi_max_np1-xi_target| = %.16g <= ||(%.16g,%.16g,%.16g)||_1\n",
		      fabs(xi_max_np1 - xi_target), p->xi[0], p->xi[1], p->xi[2]);
	      */
	      
	      load_field_variables_at_xi(p->owning_elem_id, p->xi);
	      for(i = 0; i < pdim; i++)
		current_v[i] = (1.0 - time_fraction) * (fv_old->v[i]) + time_fraction * (fv->v[i]);
	      current_Enorm = (1.0 - time_fraction) * (fv_old->Enorm) + time_fraction * (fv->Enorm);
	      fill_my_fv_old();
	      for(i = 0; i < pdim; i++)
		current_grad_Enorm[i] = (1.0 - time_fraction) * (my_fv_old->grad_Enorm[i]) + time_fraction * (fv->grad_Enorm[i]);

	      for(i = 0; i < pdim; i++)
		{
		  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
		  b[pdim+i] = mass * (p->v[i] - p->v_old[i]) - particle_dt *
		    (stokes_force_coeff * (current_v[i] - p->v[i]) +
		     dielectrophoretic_force_coeff * 2.0 * current_Enorm * current_grad_Enorm[i]);
		  /*
		  fprintf(stderr, "XI_MAX b[%d] = %g * (%g) - %g * ( %g * (%g) + %g * %g * %g\n",
			  i, mass, p->v[i]-p->v_old[i], particle_dt, stokes_force_coeff, current_v[i] - p->v[i],
			  dielectrophoretic_force_coeff*2.0, current_Enorm, current_grad_Enorm[i]);
		  */
		}
	      b_norm = nnorm(2*pdim, b);

	      newton_iterations_2++;
	    }
	  
	  if(newton_iterations_2 == 100 &&
	     b_norm > 1.0e-8 &&
	     update_norm > 1.0e-8)
	    {
	      fprintf(stderr, "\nXI_MAX move_particle() says:\n");
	      fprintf(stderr, "XI_MAX fv->x = (%g,%g,%g)\n",
		      fv->x[0], fv->x[1], fv->x[2]);
	      fprintf(stderr, "XI_MAX update = (%g,%g,%g,%g,%g,%g)\n", -x_sol[0], -x_sol[1], -x_sol[2], -x_sol[3], -x_sol[4], -x_sol[5]);
	      fprintf(stderr, "XI_MAX resid = (%g,%g,%g,%g,%g,%g), resid_norm = %g\n",
		      b[0], b[1], b[2], b[3], b[4], b[5], b_norm);
	      fprintf(stderr, "XI_MAX particle_dt = %g\n", particle_dt);
	      dump1(EXIT, p, "XI_MAX Failed to converge in 100 Newton iterations.");
	    }
	}

      return_value = newton_iterations + newton_iterations_2;
      break;

    case INERTIAL_TRACER_EXPLICIT:
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      
      for(i = 0; i < pd->Num_Dim; i++)
	{
	  total_velocity[i] = p->v[i] + particle_dt *
	    stokes_force_coeff/mass * (current_v[i] - p->v[i]);
	}
      
      if(Num_Impermeable_PBCs)
	enforce_impermeable_PBCs(p, total_velocity, particle_dt);

      for(i = 0; i < pd->Num_Dim; i++)
	{
	  p->x[i] += total_velocity[i] * particle_dt;
	  p->v[i] = total_velocity[i];
	}

      return_value = 0;
      break;

    case TRACER_IMPLICIT:
      /* Implicitly move particles at the fluid velocity. */
      memset(b, 0, 4 * sizeof(dbl));
      memset(x_sol, 0, 4 * sizeof(dbl));
      b[0] = p->x[0] - p->x_old[0] - particle_dt * current_v[0];
      b[1] = p->x[1] - p->x_old[1] - particle_dt * current_v[1];
      b_norm = nnorm(2, b);

      update_norm = 1.0e+10;
      newton_iterations = 0;
      while(b_norm > 1.0e-8 &&
	    update_norm > 1.0e-8 &&
	    newton_iterations < 100)
	{
	  memset(J_2x2, 0, 4 * sizeof(dbl));
	  /* See comments about fv->grad_v on
	   * CHARGED_TRACER_IMPLICIT. */
	  J_2x2[0][0] = 1.0 - particle_dt * fv->grad_v[0][0];
          J_2x2[0][1] = -particle_dt * fv->grad_v[1][0];
	  J_2x2[1][0] = -particle_dt * fv->grad_v[0][1];
	  J_2x2[1][1] = 1.0 - particle_dt * fv->grad_v[1][1];

	  solve_NxN_system(&J_2x2[0][0], b, x_sol, 2, 2);

	  update_norm = nnorm(2, x_sol);
	  /* include min_side_length, too? */
	  relax_factor = 1.0 - (update_norm / (update_norm + nnorm(2, current_v)));
	  if(relax_factor < 0.1 && relax_factor != 0.0)
	    fprintf(stderr, "NEWTON Relax factor = %g\n", relax_factor);
	  if(relax_factor == 0.0)
	    relax_factor = 1.0;
	  for(i = 0; i < 2; i++)
	    x_sol[i] *= relax_factor;
	  p->x[0] -= x_sol[0];
	  p->x[1] -= x_sol[1];

	  if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
	    {
	      fprintf(stderr, "\nmove_particle() says:\n");
	      fprintf(stderr, "last update = (%g,%g), resid = (%g,%g)\nresid_norm = %g\n",
		      -x_sol[0], -x_sol[1], b[0], b[1], b_norm);
	      fprintf(stderr, "particle_dt = %g\n", particle_dt);
	      dump1(EXIT, p, "Could not find element coordinates.");;
	    }
	  load_field_variables_at_xi(p->owning_elem_id, p->xi);
	  for(i = 0; i < DIM; i++)
	    current_v[i] = (1.0 - time_fraction) * (fv_old->v[i]) + time_fraction * (fv->v[i]);

	  b[0] = p->x[0] - p->x_old[0] - particle_dt * current_v[0];
	  b[1] = p->x[1] - p->x_old[1] - particle_dt * current_v[1];
	  b_norm = nnorm(2, b);

	  newton_iterations++;
	}

      if(newton_iterations == 100 &&
	 b_norm > 1.0e-8 &&
	 update_norm > 1.0e-8)
	dump1(EXIT, p, "Too many Newton iterations.");

      memcpy(p->v, current_v, DIM * sizeof(dbl));
      return_value = newton_iterations;
      break;

    case INERTIAL_TRACER_IMPLICIT:
      /* These need to match the same values used in compute_particle_dt(). */
      /* fluid_mass = 4.0/3.0 * M_PIE * mp->density * rcubed; */
      stokes_force_coeff = 6.0 * M_PIE * mp->viscosity * Particle_Radius;
      gravity_force_coeff = Particle_Model_Data[1];

      for(i = 0; i < pdim; i++)
	vel_rel[i] = fv->v[i] - p->v[i];
      vel_diff = nnorm(pdim, vel_rel);
      Re_p = Particle_Density * vel_diff * 2.0 * Particle_Radius / mp->density;
      if(Re_p <= 0.1)
	Re_p_correction = 1.0 + 3.0/16.0 * Re_p;
      else
	Re_p_correction = 1.0 + 0.0565 * pow(Re_p, 0.525);

      if(Re_p > 500.0)
	fprintf(stderr, "Applying possibly improper Re_p_correction, Re_p = %g\n", Re_p);
      stokes_force_coeff *= Re_p_correction;

      /*
      for(i = 0; i < pdim; i++)
	{
	  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	  b[pdim+i] = mass * (p->v[i] - p->v_old[i]) - particle_dt *
	    (stokes_force_coeff * (current_v[i] - p->v[i]));
	}
      b[2*pdim-1] -= particle_dt * gravity_force_coeff * (mass - fluid_mass);
      b_norm = nnorm(2*pdim, b);
      */

      /* Forces are Stokes, gravity, and pressure gradient.  The latter two are often combined (for vertical flow) as
	 "buoyancy".

	 Gravity force = (bubble volume) * (bubble density) * g, where g is a vector, e.g., -980(0,1).

	 Pressure force = - (bubble volume) * grad(p).

	 If p = - (fluid density) * g * h, then grad(p) = - (fluid density) * g.
	 Thus, pressure force = - (bubble volume) * (fluid density) * g.
	 Also, gravity force = (bubble volume) * (bubble density) * g.
	 And so if fluid density = bubble density, total force = 0.

	 Note that the interesting cases happen when grad(p) is not parallel to (0,1).
      */

      for(i = 0; i < pdim; i++)
	{
	  b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	  b[pdim+i] = p->v[i] - p->v_old[i] - particle_dt / mass *
	    (stokes_force_coeff * (current_v[i] - p->v[i]) +
	     -fv->grad_P[i] * volume);
	}
      b[gravity_index] -= particle_dt * gravity_force_coeff;

      b_norm = nnorm(2*pdim, b);

      update_norm = 1.0e+10;
      newton_iterations = 0;
      while(b_norm > 1.0e-8 &&
	    update_norm > 1.0e-8 &&
	    newton_iterations < 100)
	{
	  memset(J, 0, 36 * sizeof(dbl));
	  /* See comments about fv->grad_v on
	   * CHARGED_TRACER_IMPLICIT. */
	  for(i = 0; i < pdim; i++)
	    {
	      J[i][i] = 1.0;
	      J[i][pdim+i] = -particle_dt;
	      J[pdim+i][pdim+i] = 1.0 + particle_dt * stokes_force_coeff/mass;
	      for(j = 0; j < pdim; j++)
		J[pdim+i][j] = -particle_dt * stokes_force_coeff/mass * fv->grad_v[j][i];
	    }

	  solve_NxN_system(&J[0][0], b, x_sol, 2*pdim, 6);
	  
	  update_norm = nnorm(2*pdim, x_sol);

	  /* include min_side_length, too? */
	  relax_factor = 1.0 - (update_norm / (update_norm + nnorm(pdim, p->v)));
	  if(nnorm(pdim, p->v) <= 1.0e-8)
	    relax_factor = 1.0;
	  if(relax_factor < 0.1 && relax_factor != 0.0)
	    fprintf(stderr, "NEWTON Relax factor = %g\n", relax_factor);
	  
	  for(i = 0; i < 2*pdim; i++)
	    x_sol[i] *= relax_factor;

	  for(i = 0; i < pdim; i++)
	    {
	      p->x[i] -= x_sol[i];
	      p->v[i] -= x_sol[pdim+i];
	    }

	  if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
	    dump1(EXIT, p, "Could not find element coordinates.");

	  if(p->xi[0] < -1.5) p->xi[0] = -1.5;
	  if(p->xi[1] < -1.5) p->xi[1] = -1.5;
	  if(p->xi[2] < -1.5) p->xi[2] = -1.5;
	  if(p->xi[0] > 1.5)  p->xi[0] = 1.5;
	  if(p->xi[1] > 1.5)  p->xi[1] = 1.5;
	  if(p->xi[2] > 1.5)  p->xi[2] = 1.5;

	  load_field_variables_at_xi(p->owning_elem_id, p->xi);

	  for(i = 0; i < pdim; i++)
	    current_v[i] = (1.0 - time_fraction) * fv_old->v[i] + time_fraction * fv->v[i];

	  for(i = 0; i < pdim; i++)
	    {
	      b[i] = p->x[i] - p->x_old[i] - particle_dt * p->v[i];
	      b[pdim+i] = p->v[i] - p->v_old[i] - particle_dt / mass *
		(stokes_force_coeff * (current_v[i] - p->v[i]) +
		 -fv->grad_P[i] * volume);
	    }
	  b[gravity_index] -= particle_dt * gravity_force_coeff;

	  b_norm = nnorm(2*pdim, b);
	  
	  newton_iterations++;
	}

      if(newton_iterations == 100 &&
	 b_norm > 1.0e-8 &&
	 update_norm > 1.0e-8)
	{
	  fprintf(stderr, "\nmove_particle() says:\n");
	  fprintf(stderr, "fv->x = (%g,%g,%g)\n", fv->x[0], fv->x[1], fv->x[2]);
	  fprintf(stderr, "fv->current_v = (%g,%g,%g)\n", current_v[0], current_v[1], current_v[2]);
	  fprintf(stderr, "update = (%g,%g,%g,%g,%g,%g)\n", -x_sol[0], -x_sol[1], -x_sol[2], -x_sol[3], -x_sol[4], -x_sol[5]);
	  fprintf(stderr, "resid = (%g,%g,%g,%g,%g,%g), resid_norm = %g\n",
		  b[0], b[1], b[2], b[3], b[4], b[5], b_norm);
	  fprintf(stderr, "p->v[i] -  p->v_old[i] - particle_dt/mass * (stokes * (current_v[i] - p->v[i]) - fv->grad_P[i] * volume) =\n");
	  fprintf(stderr, "  i=0: %g -  %g - %g / %g * (%g * (%g - %g) - %g * %g)\n", p->v[0], p->v_old[0], particle_dt, mass, stokes_force_coeff, current_v[0], p->v[0], fv->grad_P[0], volume);
	  fprintf(stderr, "  i=1: %g -  %g - %g / %g * (%g * (%g - %g) - %g * %g)\n", p->v[1], p->v_old[1], particle_dt, mass, stokes_force_coeff, current_v[1], p->v[1], fv->grad_P[1], volume);
	  fprintf(stderr, "  i=2: %g -  %g - %g / %g * (%g * (%g - %g) - %g * %g)\n", p->v[2], p->v_old[2], particle_dt, mass, stokes_force_coeff, current_v[2], p->v[2], fv->grad_P[2], volume);
	  fprintf(stderr, "b[%d] -= particle_dt * gravity =\n", gravity_index);
	  fprintf(stderr, "  b[%d] -= %g * %g\n", gravity_index, particle_dt, gravity_force_coeff);
	  fprintf(stderr, "fv->v = (%g,%g,%g), fv_old->v = (%g,%g,%g)\n",
		  fv->v[0], fv->v[1], fv->v[2], fv_old->v[0], fv_old->v[1], fv_old->v[2]);
	  fprintf(stderr, "particle_dt = %g\n", particle_dt);
	  dump1(EXIT, p, "Failed to converge in 100 Newton iterations.");
	}
      return_value = newton_iterations;
      break;

    default:
      dump1(EXIT, p, "Unknown Particle_Model.");
    }

  if(Particle_Move_Domain == BRICK)
    for(i = 0; i < pdim; i++)
      {
	p->x[i] = MAX(p->x[i], Particle_Move_Domain_Reals[2*i]);
	p->x[i] = MIN(p->x[i], Particle_Move_Domain_Reals[2*i+1]);
      }

  *max_particle_dt = particle_dt;
  return return_value;
}


/* Return distance between two points, assuming dimension of points is
 * my_dim. */
static dbl
my_distance(const dbl *x,
	    const dbl *y,
	    int my_dim)
{
  int i;
  dbl tmp, tmp2;

  tmp = 0.0;
  for(i = 0; i < my_dim; i++)
    {
      tmp2 = x[i] - y[i];
      tmp += tmp2 * tmp2;
    }
  return sqrt(tmp);
}


/* This routine modifies the particles' pre-existing total_velocity
 * vector if the particle is within the appropriate distance from an
 * impermeable particle boundary condition. */
static void
enforce_impermeable_PBCs(particle_t * const p,
			 dbl * total_velocity,
			 dbl particle_dt)
{
}


/* This routine returns a random number from a standard normal
 * distribution.  Since the basic algoritm actually computes two
 * random numbers each time, we use some static data to save us half
 * the work.  From Numerical Recipes (ugh).
 *
 * This is *not* a truncated normal distribution. */
static dbl
drand_standard_normal(void)
{
  static int have_spare = 0;
  static dbl x_save = 0.0;
  dbl pt[2], pt_norm, factor;
  pt[0] = 0;
  pt[1] = 0;

  if(have_spare)
    {
      have_spare = 0;
      return x_save;
    }
  else
    {
      pt_norm = 0.0;
      while(pt_norm >= 1.0 || pt_norm == 0.0)
	{
	  pt[0] = 2.0 * drand48() - 1.0;
	  pt[1] = 2.0 * drand48() - 1.0;
	  pt_norm = pt[0]*pt[0] + pt[1]*pt[1];
	}
      factor = sqrt(-2.0 * log(pt_norm) / pt_norm);
      x_save = pt[1] * factor;
      have_spare = 1;
      return pt[0] * factor;
    }
}


/* This is a truncated normal distribution where the result is
 * restricted to be within +/- 3 standard deviations.  The > 3*sigma
 * results should not be significant in the overall numerical coupling
 * (or you're probably doing something very unstable/stiff/etc.)  So
 * our "stochastic noise" won't cause particles to (eventually) exceed
 * the speed of light, etc... */
static dbl
drand_truncated_normal(const dbl mean,
		       const dbl stddev)
{
  dbl x = 3.1;

  while(fabs(x) > 3.0)
    x = drand_standard_normal();
  return mean + stddev * x;
}


/* This is the main dump and-possibly-exit routine.  It spits out all
 * the particle information we have, as well as the current static
 * values (for element coordinates loaded, etc.).  It allows the
 * calling code to include special information in the variable
 * argument list that varies from call to call.
 *
 * Unfortunately, there is no way to do a count of the arguments in
 * the variable argument list, so the minimum call looks like
 * dump(EXIT, p, ""); -- that last "" cannot be left off.
 *
 * Furthermore, to make this compatible with the macro definition and
 * the problem with including/excluding the variable argument list,
 * the first "variable arugment" will always be the format line...
 */
static void
dump_fcn(const char * const file_name,
	 const int line_number,
	 const int exit_too,
	 particle_t * p,
	 ...)
{
  va_list args;
  char *fmt, state_string[80];
  int i;

  switch(p->state)
    {
    case ACTIVE: strcpy(state_string, "ACTIVE"); break;
    case DEAD: strcpy(state_string, "DEAD"); break;
    case PROC_TRANSFER: strcpy(state_string, "PROC_TRANSFER"); break;
    case TERMINATION_MARKER: strcpy(state_string, "TERMINATION_MARKER"); break;
    case GHOST: strcpy(state_string, "GHOST"); break;
    default: strcpy(state_string, "<UNKNOWN>"); break;
    }

  fprintf(stderr, "\n=== DUMP BEGINS ===\n");
#ifdef PARALLEL
  fprintf(stderr, "Called from %s:%d on Proc %d\n", file_name, line_number, ProcID);
#else
  fprintf(stderr, "Called from %s:%d\n", file_name, line_number);
#endif

  fprintf(stderr, "\nParticle information\n--------------------\n");
  fprintf(stderr, "  p->x = (%g,%g,%g), p->x_old = (%g,%g,%g)\n",
	  p->x[0], p->x[1], p->x[2], p->x_old[0], p->x_old[1], p->x_old[2]);
  fprintf(stderr, "  p->xi = (%g,%g,%g), p->xi_old = (%g,%g,%g)\n",
	  p->xi[0], p->xi[1], p->xi[2], p->xi_old[0], p->xi_old[1],p->xi_old[2]);
  fprintf(stderr, "  p->v = (%g,%g,%g), p->v_old = (%g,%g,%g)\n",
	  p->v[0],p->v[1], p->v[2], p->v_old[0], p->v_old[1], p->v_old[2]);
  fprintf(stderr, "  p->time = %g, p->time_old = %g\n",
	  p->time, p->time_old);
  fprintf(stderr, "  p->state = %s, p->owning_elem_id = %d, p->output_sample_number = %d\n",
	  state_string, p->owning_elem_id, p->output_sample_number);

  fprintf(stderr, "\nStatic information\n------------------\n");
  fprintf(stderr, "  last_element_loaded = %d\n",
	  last_element_loaded);
  fprintf(stderr, "  last_elements_nodes_loaded = %d\n",
	  last_elements_nodes_loaded);
  fprintf(stderr, "  node coordinates =\n");
  for(i = 0; i < nodes_per_element; i++)
    fprintf(stderr, "    [%d] = (%g,%g,%g)\n",
	    i, node_coord[i][0], node_coord[i][1], node_coord[i][2]);
  fprintf(stderr, "  last_elements_fv_loaded = %d\n",
	  last_elements_fv_loaded);
  fprintf(stderr, "  last_fvs_xi_loaded = (%g,%g,%g)\n",
	  last_fvs_xi_loaded[0], last_fvs_xi_loaded[1], last_fvs_xi_loaded[2]);

  va_start(args, p);
  fmt = (char *) va_arg(args, char *);
  if(strlen(fmt))
    {
      fprintf(stderr, "\nAdditional information\n----------------------\n");
      vfprintf(stderr, fmt, args);
      fprintf(stderr, "\n");
    }
  va_end(args);
  fflush(stderr);

  if(exit_too)
    {
#ifdef PARALLEL
      MPI_Finalize();
#endif
      exit(-1);
    }
}


/* I need some unknown values at the old timestep.  Some of these
 * aren't available in the fv_old structure.  Instead of adding them
 * to the Diet_Field_Variables structure (which is on a diet), I'll
 * just cook up local copies here.  NOTE: these derivatives are *not*
 * correct for a deforming mesh.  They would need to be corrected with
 * the xdot value, re: fv_old->v = fv->v - xdot . fv->grad_v.  This
 * assumes the proper element has already been loaded.
 */
static void
fill_my_fv_old(void)
{
  int i, j, dofs;
  
  if(pd->v[pg->imtrx][VOLTAGE])
    {
      dofs  = ei[pg->imtrx]->dof[VOLTAGE];
      for(i = 0; i < VIM; i++)
	{
	  my_fv_old->grad_V[i] = 0.0;
	  for(j = 0; j < dofs; j++)
	    my_fv_old->grad_V[i] += *esp_old->V[j] * bf[VOLTAGE]->grad_phi[j][i];
	}
    }

  if(pd->v[pg->imtrx][ENORM])
    {
      dofs  = ei[pg->imtrx]->dof[ENORM];
      for(i = 0; i < VIM; i++)
	{
	  my_fv_old->grad_Enorm[i] = 0.0;
	  for(j = 0; j < dofs; j++)
	    my_fv_old->grad_Enorm[i] += *esp_old->Enorm[j] * bf[ENORM]->grad_phi[j][i];
	}
    }
}


/* This is the main driver for advancing a particle from time time
 * step N to goma time step N+1.  It and its children increment the
 * time accumulators (currently output, particle (move), and
 * communication).  It also increments the max_newton_iterations and
 * max_particle_iterations if necessary.  It is possible for
 * find_exit_wound to remove this particle from consideration, so !p
 * may become true "prematurely".
 *
 * I've set it up so that the caller is responsible for add/removing
 * particles from elements lists and processor lists...  Don't do that
 * here or things will break.
 */
static void
advance_a_particle(particle_t * p,
		   const dbl global_start_time,
		   const dbl global_end_time,
		   const int n)
{
  dbl particle_dt, output_start_ust;
  int newton_iterations, particle_iterations;
  /*  int original_owning_elem_id; */
  int done;
  
  load_element_information(p->owning_elem_id);

  max_newton_iterations = 0;
  particle_iterations = 0;

  done = 0;

  if(p->time >= global_end_time)
    {
      /* I probably got here b/c I was a new particle that finished
       * its last timestep in the transfer work... */
      done = 1;
    }

  while(!done)
    {
      memcpy(p->v_old, p->v, DIM * sizeof(dbl));
      memcpy(p->x_old, p->x, DIM * sizeof(dbl));
      memcpy(p->real_data_old, p->real_data, MAX_DATA_REAL_VALUES * sizeof(dbl));
      memcpy(p->int_data_old, p->int_data, MAX_DATA_INT_VALUES * sizeof(int));

      /* What is the velocity (and other unknowns) at our particle's
       * location?  We computed xi at the end of the last timestep. */
      load_field_variables_at_xi(p->owning_elem_id, p->xi);

      /* Get the time-step size. */
      p->time_old = p->time;
      particle_dt = compute_particle_dt(p, global_end_time - p->time);

      /* Move the particle.  Lots of different methods under the
       * hood here, depending on the Particle_Model. */
      newton_iterations = move_particle(p, &particle_dt, global_start_time, global_end_time - global_start_time);
      p->time += particle_dt;

      /*
      fprintf(stderr, "NEW XI = (%g,%g,%g)\n", p->xi[0], p->xi[1], p->xi[2]);
      */
      
      max_newton_iterations = MAX(max_newton_iterations, newton_iterations);

      /* Update the element coordinates and owning elements for this particle */
      memcpy(p->xi_old, p->xi, DIM * sizeof(dbl));
	  
      if(get_element_xi_newton(p->owning_elem_id, p->x, p->xi) == -1)
	dump2(EXIT, p, "particle_dt = %g\nCould not get particle coordinates.", particle_dt);

      /* Find the particles that have moved out of their element.
       * Some of these have moved completely out of the computational
       * domain. */
      /* original_owning_elem_id = p->owning_elem_id; */
      if(fabs(p->xi[0]) >= 1.0 || fabs(p->xi[1]) >= 1.0 || fabs(p->xi[2]) >= 1.0)
	{
#ifdef SHOW_PARTICLE_MOVEMENT
	  fprintf(stderr, "PA");
#endif
	  /*
	  fprintf(stderr, "Proc%d: Calling few() from 2.\n", ProcID);
	  */
	  find_exit_wound(p->owning_elem_id, p, p->x_old, p->x, p->xi, 0, 0);
	  /*
	  fprintf(stderr, "Proc%d: that particle is state %d\n", ProcID, p->state);
	  */
	}

      /* Save the requested data for output later. */
      if(p->output_sample_number != -1)
	store_particle_data(p);
      
      /*
      fprintf(stderr, "Proc%d: outputting for particle at (%g,%g), p->time = %g, p->time_old = %g\n",
	      ProcID, p->x[0], p->x[1], p->time, p->time_old);
      */
      
      output_start_ust = ust();
      if(TimeIntegration == STEADY)
	{
	  if(p->time == global_end_time) /* because end_time is an output time in this case. */
	    output_a_particle(p, p->time_old, p->time, n, 1, 0);
	}
      else
	output_a_particle(p, p->time_old, p->time, n, (int)(p->time == global_end_time), 0);
      output_accum_ust += MAX(ust() - output_start_ust, 0.0);

      /*
      if(p->state == ACTIVE &&
	 p->owning_elem_id == original_owning_elem_id)
	{
	  if(p->output_sample_number != -1)
	    store_particle_data(p);
	  
	  output_start_ust = ust();
	  if(TimeIntegration == STEADY)
	    {
	      if(p->time == global_end_time)
		output_a_particle(p, p->time_old, p->time, n, 1, 0);
	    }
	  else
	    output_a_particle(p, p->time_old, p->time, n, (int)(p->time == global_end_time), 0);
	  output_accum_ust += max(ust() - output_start_ust, 0.0);
	}
      else
	done = 1;
      */
      
      if(p->time == global_end_time ||
	 p->state != ACTIVE)
	done = 1;

      particle_iterations++;
    }

  max_particle_iterations = MAX(max_particle_iterations, particle_iterations);
}


/* This routine removes a particle from the element_particle_list_head
 * headed list.  It fixes the pointers to keep things consistent. */
static void
remove_from_element_particle_list(particle_t * p,
				  const int elem_id)
{
  if(!(p->last))	/* If I'm the head */
    element_particle_list_head[elem_id] = p->next;
  else		/* Fix my previous neighbor */
    p->last->next = p->next;
  
  if(p->next)	/* Fix my next neighbor */
    p->next->last = p->last;
}
				 

/* This routine moves a particle to the outgoing send list.  Once it
 * is actually sent its space will be free()'ed. */
static void
add_to_send_list(particle_t * p)
{

  if(!particles_to_send)
    {
      particles_to_send = p;
      p->last = NULL;
      p->next = NULL;
    }
  else
    {
      particles_to_send->last = p;
      p->next = particles_to_send;
      particles_to_send = p;
      p->last = NULL;
    }
}


/* This routine moves a particle to the incoming list.  Space needs to
 * be allocated, too. */
static void
add_to_do_list(particle_t * p)
{
  particle_t * p_recv;

  if(!particles_to_do)
    {
      particles_to_do = (particle_t *)malloc(sizeof(particle_t));
      if(!particles_to_do)
	EH(-1, "Could not malloc particle space.");
      p_recv = particles_to_do;
      memcpy(p_recv, p, sizeof(particle_t));
      p_recv->last = NULL;
      p_recv->next = NULL;
    }
  else
    {
      p_recv = (particle_t *)malloc(sizeof(particle_t));
      if(!p_recv)
	EH(-1, "Could not malloc particle space.");
      memcpy(p_recv, p, sizeof(particle_t));
      particles_to_do->last = p_recv;
      p_recv->next = particles_to_do;
      particles_to_do = p_recv;
      p_recv->last = NULL;
    }
  p_recv->state = ACTIVE;
}


/* This routine handles whatever needs to be saved off, etc., to
 * influence the continuum solution with repsect to the particles'
 * presence.  It assumes that all particles contribute (irresepective
 * of possible p->state settings). */
static void
couple_to_continuum(void)
{
  dbl source_mass = 0.0, total_particles = 0.0;
  int source_eqn = -1;
  int i, j;
  particle_t *p;

#ifdef PARALLEL
  double temp_ust, local_total_particles;
  int recipient_proc, recipient_local_elem_id, local_num_particle_ghosts, use_ghost_particles;
  int *num_particles_in_elem, *num_ghost_particles_to_send, *num_ghost_particles_to_recv, *num_ghost_particles_copied;
  particle_t *p_tmp, **ghost_particles_to_send, *ghost_particles_to_recv;
  MPI_Status mpi_status;

  use_ghost_particles = (Particle_Model == SWIMMER_EXPLICIT ||
			 Particle_Model == SWIMMER_IMPLICIT);

  /* Testing to see if we really do need to ghost particles, or if
   * they already get summed properly in the source terms now,
   * anyways.  After thinking about it I think that they might
   * actually get summed properly, already!. */
  use_ghost_particles = 1;

  local_num_particle_ghosts = 0;

  if(use_ghost_particles)
    {
      /* Exchange ghost particles first! */
      num_ghost_particles_to_send = (int *)calloc((unsigned)Num_Proc, sizeof(int));
      num_ghost_particles_to_recv = (int *)calloc((unsigned)Num_Proc, sizeof(int));
      ghost_particles_to_send = (particle_t **)calloc((unsigned)Num_Proc, sizeof(particle_t *));
      num_particles_in_elem = (int *)calloc((unsigned)(static_exo->num_elems), sizeof(int));
      num_ghost_particles_copied = (int *)calloc((unsigned)Num_Proc, sizeof(int));
 
      /* I should introduce a level of indirection into
       * element_particle_info[i].ghost_proc[] if I want to speed things
       * up.  I think it doesn't matter, really... */

      /* First count the number of particles currently in each of the
       * to-be-ghosted elements. */
      for(i = 0; i < static_exo->num_elems; i++)
	if(element_particle_info[i].num_ghost_target_elems)
	  for(p = element_particle_list_head[i]; p; p = p->next)
	    num_particles_in_elem[i]++;

      /* Now count the number of particles going to each of the other
       * processors. */
      for(i = 0; i < static_exo->num_elems; i++)
	for(j = 0; j < element_particle_info[i].num_ghost_target_elems; j++)
	  num_ghost_particles_to_send[element_particle_info[i].ghost_proc[j]] += num_particles_in_elem[i];
      for(i = 0; i < Num_Proc; i++)
	local_num_particle_ghosts += num_ghost_particles_to_send[i];

      /* Allocate the space for the outgoing particles. */
      for(i = 0; i < Num_Proc; i++)
	if(num_ghost_particles_to_send[i])
	  ghost_particles_to_send[i] = (particle_t *)calloc((unsigned)num_ghost_particles_to_send[i], sizeof(particle_t));
	else
	  ghost_particles_to_send[i] = NULL;
      
      /* Now memcpy the particle info into the send buffer.  Set the
       * owning element id, too! */
      for(i = 0; i < static_exo->num_elems; i++)
	if(element_particle_info[i].num_ghost_target_elems)
	  for(p = element_particle_list_head[i]; p; p = p->next)
	    for(j = 0; j < element_particle_info[i].num_ghost_target_elems; j++)
	      {
		recipient_proc = element_particle_info[i].ghost_proc[j];
		recipient_local_elem_id = element_particle_info[i].ghost_local_elem_id[j];
		memcpy(&ghost_particles_to_send[recipient_proc][num_ghost_particles_copied[recipient_proc]], p, sizeof(particle_t));
		ghost_particles_to_send[recipient_proc][num_ghost_particles_copied[recipient_proc]].owning_elem_id = recipient_local_elem_id;
		num_ghost_particles_copied[recipient_proc]++;
	      }

      /* Now we round-robin the information about... */
      for(i = 0; i < Num_Proc; i++)
	{
	  if(i == ProcID)
	    {
	      temp_ust = ust();
	      MPI_Bcast(num_ghost_particles_to_send, Num_Proc, MPI_INT, i, MPI_COMM_WORLD);
	      for(j = 0; j < Num_Proc; j++)
		if(i != j && num_ghost_particles_to_send[j])
		  MPI_Send(&ghost_particles_to_send[j][0], num_ghost_particles_to_send[j] * sizeof(particle_t), MPI_BYTE, j, 0, MPI_COMM_WORLD);
	
	      communication_accum_ust += MAX(ust() - temp_ust, 0.0);
	    }
	  else
	    {
	      temp_ust = ust();
	      MPI_Bcast(num_ghost_particles_to_recv, Num_Proc, MPI_INT, i, MPI_COMM_WORLD);
	      communication_accum_ust += MAX(ust() - temp_ust, 0.0);
	      if(num_ghost_particles_to_recv[ProcID])
		{
		  ghost_particles_to_recv = (particle_t *)calloc((unsigned)num_ghost_particles_to_recv[ProcID], sizeof(particle_t));
		  temp_ust = ust();
		  MPI_Recv(ghost_particles_to_recv, num_ghost_particles_to_recv[ProcID] * sizeof(particle_t), MPI_BYTE, i, 0, MPI_COMM_WORLD, &mpi_status);
		  communication_accum_ust += MAX(ust() - temp_ust, 0.0);
		  for(j = 0; j < num_ghost_particles_to_recv[ProcID]; j++)
		    {
		      ghost_particles_to_recv[j].state = GHOST;
		      p = create_a_particle(&ghost_particles_to_recv[j], ghost_particles_to_recv[j].owning_elem_id);
		      num_particles--; /* Don't count these as new particles, they are counted as ghosts. */
		    }
		  free(ghost_particles_to_recv);
		}
	    }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }
  
  total_num_particle_ghosts = 0;
  MPI_Reduce(&local_num_particle_ghosts, &total_num_particle_ghosts, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  
  /* Model-specific initializations. */
  if(Particle_Model == SWIMMER_EXPLICIT ||
     Particle_Model == SWIMMER_IMPLICIT)
    {
      source_eqn = R_MOMENTUM1 + pdim - 1;
      if(mp->DensityModel != CONSTANT)
	EH(-1, "Sorry, only handling CONSTANT density models right now.  Check back later.");
      source_mass = 4.0/3.0 * M_PIE  * (Particle_Density - mp->density) * Particle_Radius * Particle_Radius * Particle_Radius * Particle_Ratio;
    }

  for(i = 0; i < static_exo->num_elems; i++)
    {
      for(p = element_particle_list_head[i]; p; p = p->next)
	{
	  if(Particle_Model == SWIMMER_EXPLICIT ||
	     Particle_Model == SWIMMER_IMPLICIT)
	    {
	      load_field_variables_at_xi(i, p->xi);
		  
	      if(p->owning_elem_id != i)
		{
		  fprintf(stderr, "PROC%d: UH-OH, THEY ARE NOT EQUAL!\n", ProcID);
		  fprintf(stderr, "PROC%d: \t\ti = %d, particle->owning_elem_id = %d\n", ProcID, i, p->owning_elem_id);
		  fprintf(stderr, "PROC%d: \t\tstate = %d\n", ProcID, p->state);
		}

	      for(j = 0; j < ei[pg->imtrx]->dof[source_eqn]; j++)
		{
		  element_particle_info[p->owning_elem_id].source_term[j] += source_mass * bf[source_eqn]->phi[j];
		  total_particles += bf[source_eqn]->phi[j];
		}
	    }
	}
    }

#ifdef PARALLEL
  local_total_particles = total_particles;
  MPI_Reduce(&local_total_particles, &total_particles, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

#ifdef PARALLEL
  /*
  if(ProcID == 0)
    fprintf(stderr, "    (TOTAL PARTICLES (+GHOSTS) = %g)\n", total_particles);
  */
    
  if(use_ghost_particles)
    {
      /* Currently there is no need to keep the ghost particles around,
       * but there might be one day... */
      for(i = 0; i < static_exo->num_elems; i++)
	{
	  p = element_particle_list_head[i];
	  while(p)
	    {
	      p_tmp = p->next;
	      if(p->state == GHOST)
		{
		  remove_from_element_particle_list(p, i);
		  free(p);
		}
	      p = p_tmp;
	    }
	}
      for(i = 0; i < Num_Proc; i++)
	if(ghost_particles_to_send[i])
	  free(ghost_particles_to_send[i]);
      free(num_ghost_particles_to_send);
      free(num_ghost_particles_to_recv);
      free(num_ghost_particles_copied);
      free(ghost_particles_to_send);
    }
#endif
}

/*
 * This routine is used to test the integrity of the element<->element
 * map.  It pretends to move a particle from the centroid of the
 * element through the centriod of each of the 6 sides.  This assumes
 * that (1) we're in 3D and (2) that the elements are not too
 * distorted so the set of exit vectors are useful.
 *
 * I'm an idiot.  I don't really need all that vector stuff, but I'm
 * leaving it in here in case it is useful later...
 */
static void
test_map_integrity(void)
{
  int elem_id, neighbor_elem_id, neighbor_elem_proc_id, face_id, PBC_id, proc_id;
  FILE *fp;

#ifdef PARALLEL
  for(proc_id = 0; proc_id < Num_Proc; proc_id++)
    {
      if(proc_id == ProcID)
	{
	  if(proc_id == 0)
	    fp = fopen("map.txt", "w");
	  else
	    fp = fopen("map.txt", "a");
	  for(elem_id = 0; elem_id < static_exo->num_elems; elem_id++)
	    {
	      if(DPI_ptr->elem_owner[elem_id] != ProcID)
		{
		  fprintf(fp, "Proc %2d Elem %4d(%3d) <not owned>\n", ProcID, DPI_ptr->elem_index_global[elem_id], elem_id);
		  continue;
		}
	      else
		fprintf(fp, "Proc %2d Elem %4d(%3d)\n", ProcID, DPI_ptr->elem_index_global[elem_id], elem_id);
	      for(face_id = 0; face_id < sides_per_element; face_id++)
		{
		  PBC_id = element_particle_info[elem_id].PBC_side_id[face_id];
		  fprintf(fp, "\t\tFace %d PBC=%2d ", face_id, PBC_id);
		  if(PBC_id >= 0)
		    fprintf(fp, "boundary condition                 OK\n");
		  else
		    {
		      if(PBC_id == -1)
			{
			  neighbor_elem_id = static_exo->elem_elem_list[static_exo->elem_elem_pntr[elem_id] + face_id];
			  if(neighbor_elem_id == -1)
			    fprintf(fp, "missing local neighbor            FAIL\n");
			  else
			    fprintf(fp, "local neighbor=%4d(%3d) (%5s)   OK\n", DPI_ptr->elem_index_global[neighbor_elem_id],
				    neighbor_elem_id, (DPI_ptr->elem_owner[neighbor_elem_id] == ProcID ? "owned" : "ghost"));
			}
		      else
			{
			  neighbor_elem_proc_id = -(PBC_id+2);
			  neighbor_elem_id = element_particle_info[elem_id].owner_local_element_id[face_id];
			  fprintf(fp, "remote proc=%d neighbor=%3d        ????\n", neighbor_elem_proc_id, neighbor_elem_id);
			}
		    }
		}
	    }
	  fclose(fp);
	}
      fflush(stderr);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#else
  fprintf(stderr, "No integrity may for serial.");
  exit(-1);
#endif
}

/* This routine loads particles from the specified input file.  It
 * should be in a format that is identical to the "Full output stride"
 * specified file.  Note that *all* particles are loaded from this
 * file, so if your output file contains multiple output times you
 * will want to used a modified version. */
static void
load_restart_file(void)
{
  int i;
  particle_t p, *p_ptr;
  char garbage[255];
  FILE *fp;

  if(!(fp = fopen(construct_filename(Particle_Restart_Filename), "r")))
    {
      fprintf(stderr, "trying to open %s.\n", construct_filename(Particle_Restart_Filename));
      EH(-1, "Could not open restart file");
    }
  zero_a_particle(&p);
  while(!feof(fp))
    {
      char *fgetserr;
      for(i = 0; i < pdim; i++)
	if(fscanf(fp, "%lf ", &p.x[i]) != 1)
	  continue;		/* probable EOF */
      if(fscanf(fp, "%lf %d", &p.time, &p.owning_elem_id) != 2)
	continue;		/* probable EOF */
      fgetserr = fgets(garbage, 254, fp);
      if (fgetserr == NULL) {
	EH(-1, "Error reading line in particle restart file");
      }
      if(get_element_xi_newton(p.owning_elem_id, p.x, p.xi) == -1)
	EH(-1, "Could not place particle from restart.");
      load_field_variables_at_xi(p.owning_elem_id, p.xi);
      for(i = 0; i < pdim; i++)
	if(fabs(fv->x[i] - p.x[i]) > 1.0e-6)
	  EH(-1, "Did not match particle location on restart.");
      p_ptr = create_a_particle(&p, p.owning_elem_id);
      switch(Particle_Model)
	{
	case INERTIAL_TRACER_EXPLICIT:
	case INERTIAL_TRACER_IMPLICIT:
	case TRACER_EXPLICIT:
	case TRACER_IMPLICIT:
	case SWIMMER_IMPLICIT:
	case SWIMMER_EXPLICIT:
	case CHARGED_TRACER_EXPLICIT:
	case CHARGED_TRACER_IMPLICIT:
	case DIELECTROPHORETIC_TRACER_IMPLICIT:
	  memcpy(p_ptr->v, fv->v, DIM * sizeof(dbl));
	  break;
	default:
	  dump1(EXIT, p_ptr, "Unknown Particle_Model.  You shouldn't be here.");
	}
      if(p_ptr->output_sample_number > -1)
	store_particle_data(p_ptr);
      output_a_particle(p_ptr, 0.0, 0.0, 0, 1, 1);
    }
}
