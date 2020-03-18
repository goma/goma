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
// 
/*
 * $Id: ac_particles.h,v 5.1 2007-09-18 18:53:40 prschun Exp $
 */

/*
 * Copyright (C) 2001 Sandia National Laboratories
 */

#ifndef GOMA_AC_PARTICLES_H
#define GOMA_AC_PARTICLES_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_AC_PARTICLES_C
#define EXTERN /* do nothing */
#endif

#include <stdio.h>

#include "std.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "stdio.h"
#include "ac_stability_util.h"

#ifndef GOMA_AC_PARTICLES_C
#define EXTERN extern
#endif

#define MAX_DATA_REAL_VALUES                 3
#define MAX_DATA_INT_VALUES                  1
#define MAX_PBC_REAL_VALUES                  5
#define MAX_PBC_INT_VALUES                   5
#define MAX_DOMAIN_REAL_VALUES               6
#define MAX_PARTICLE_MODEL_DATA_VALUES       7
#define MAX_PBC_STRING_DATA_LENGTH          80
#define MAX_PARTICLE_OUTPUT_VARIABLE_LENGTH 12
#define MAX_PARTICLE_FILENAME_LENGTH        80
#define MAX_PARTICLE_STRING_LENGTH         255

#define XI_BOUNDARY_TOLERANCE0 1.0e-8
#define XI_BOUNDARY_TOLERANCE1 1.0e-12
#define XI_BOUNDARY_TOLERANCE2 1.0e-14

#define MAX_RANK 6 /* Used in the solve_NxN_system routine to avoid lots of malloc()'s. */

#define EXIT    1 /* For calling the dump() routine. */
#define NO_EXIT 0

/* This is the "right" way to do it, but the current Sun compilers
 * can't handle it, so I have all the more explicit call definisions
 * below. */
/*
#define dump(ec, p, ...) dump_fcn(__FILE__, __LINE__, ec, p, __VA_ARGS__)
*/
#define dump1(ec, p, x1) dump_fcn(__FILE__, __LINE__, ec, p, x1)
#define dump2(ec, p, x1,x2) dump_fcn(__FILE__, __LINE__, ec, p, x1,x2)
#define dump3(ec, p, x1,x2,x3) dump_fcn(__FILE__, __LINE__, ec, p, x1,x2,x3)
#define dump4(ec, p, x1,x2,x3,x4) dump_fcn(__FILE__, __LINE__, ec, p, x1,x2,x3,x4)
#define dump10(ec, p, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) dump_fcn(__FILE__, __LINE__, ec, p, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
#define dump20(ec, p, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20) dump_fcn(__FILE__, __LINE__, ec, p, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)

enum Particle_Model_t { NO_PARTICLE_MODEL,
			TRACER_EXPLICIT, TRACER_IMPLICIT,
			INERTIAL_TRACER_EXPLICIT, INERTIAL_TRACER_IMPLICIT,
			SWIMMER_EXPLICIT, SWIMMER_IMPLICIT,
                        CHARGED_TRACER_EXPLICIT, CHARGED_TRACER_IMPLICIT,
			DIELECTROPHORETIC_TRACER_IMPLICIT};

enum Particle_Output_Format_t { FLAT_TEXT, TECPLOT };

enum PBC { PBC_OUTFLOW, PBC_SOURCE, PBC_TARGET, PBC_FREESTREAM_SOURCE, PBC_IMPERMEABLE };

enum Particle_State_t { ACTIVE, DEAD, PROC_TRANSFER, TERMINATION_MARKER, GHOST };

enum Particle_Domain_t { UNRESTRICTED, BRICK, ACIS_OBJECT };

typedef struct {
  enum PBC type;
  int SS_id;
  dbl real_data[MAX_PBC_REAL_VALUES];
  int int_data[MAX_PBC_INT_VALUES];
  char string_data[MAX_PBC_STRING_DATA_LENGTH];
} PBC_t;

typedef struct {
  dbl grad_V[DIM];
  dbl grad_Enorm[DIM];
} my_fv_old_t;

typedef struct p_t {
  dbl x[DIM];
  dbl x_old[DIM];
  dbl xi[DIM];
  dbl xi_old[DIM];
  dbl v[DIM];
  dbl v_old[DIM];
  dbl theta, theta_old;
  dbl phi, phi_old;
  dbl real_data[MAX_DATA_REAL_VALUES];
  dbl real_data_old[MAX_DATA_REAL_VALUES];
  dbl time, time_old;
  enum Particle_State_t state;
  int int_data[MAX_DATA_INT_VALUES];
  int int_data_old[MAX_DATA_INT_VALUES];
  int owning_elem_id;
  int output_sample_number;
  struct p_t *next, *last;
#ifdef PARALLEL
  dbl x_start[DIM];
  dbl x_end[DIM];
  int owning_proc_id;
#endif
} particle_t;

typedef struct {
  int *PBC_side_id;
#ifdef PARALLEL
  int *owner_local_element_id;  /* On other side of face on processor boundary, what's the owner's local element id? */
  int *ghost_proc;		/* What processors and their local element id's do these particles need to be sent to? */
  int *ghost_local_elem_id;
  int num_ghost_target_elems;
#endif
  dbl *source_term;
  int *list_PBC_IMPERMEABLE;
} element_particle_info_t;

typedef char particle_variable_s[MAX_PARTICLE_OUTPUT_VARIABLE_LENGTH];
typedef char particle_filename_s[MAX_PARTICLE_FILENAME_LENGTH];
typedef char particle_s[MAX_PARTICLE_STRING_LENGTH];

extern element_particle_info_t *element_particle_info;

extern int Particle_Dynamics;
extern enum Particle_Model_t Particle_Model;
extern dbl Particle_Model_Data[MAX_PARTICLE_MODEL_DATA_VALUES];
extern int Particle_Number;
extern particle_filename_s Particle_Restart_Filename;
extern int Particle_Output_Stride;
extern dbl Particle_Output_Time_Step;
extern int Particle_Max_Time_Steps;
extern enum Particle_Output_Format_t Particle_Output_Format;
extern int Particle_Full_Output_Stride;
extern particle_filename_s Particle_Full_Output_Filename;
extern int Particle_Number_Sample_Types;
extern int *Particle_Number_Samples;
extern int *Particle_Number_Samples_Existing;
extern int *Particle_Number_Output_Variables;
extern particle_variable_s **Particle_Output_Variables;
extern particle_filename_s *Particle_Filename_Template;

extern dbl Particle_Density;
extern dbl Particle_Radius;
extern dbl Particle_Ratio;
extern int Particle_Show_Debug_Info;
extern enum Particle_Domain_t Particle_Creation_Domain;
extern enum Particle_Domain_t Particle_Move_Domain;
extern particle_filename_s Particle_Creation_Domain_Filename;
extern particle_s Particle_Creation_Domain_Name;
extern particle_filename_s Particle_Move_Domain_Filename;
extern particle_s Particle_Move_Domain_Name;
extern dbl Particle_Creation_Domain_Reals[MAX_DOMAIN_REAL_VALUES];
extern dbl Particle_Move_Domain_Reals[MAX_DOMAIN_REAL_VALUES];

extern dbl xi_boundary_tolerances[3];


extern int Particle_Number_PBCs;
extern PBC_t *PBCs;

EXTERN void initialize_particles
(const Exo_DB *,
       dbl * const,
       dbl * const,
       dbl * const,
       dbl * const,
       dbl * const);

EXTERN int compute_particles
(const Exo_DB *,
       dbl * const,
       dbl * const,
       dbl * const,
       dbl * const,
       dbl * const,
       const dbl,
       const dbl,
       const int);

EXTERN void rd_particle_specs	/* mm_input_particles.c */
(FILE *,
       char *);
#endif /* GOMA_AC_PARTICLES_H */

