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
 

#include <assert.h>
#include "mm_mp_const.h"

#define	NA "NA"

#define TIME_ORD 1
#define X_ORD    2
#define Y_ORD    3
#define YZ_ORD   4
#define YX_ORD   5
#define ZX_ORD   6
#define ZY_ORD   7
#define Z_ORD	 8
#define XY_ORD	 9
#define XZ_ORD	10

#define VELOCITY1_ABS 1
#define U_ABS 2
#define VELOCITY2_ABS 3
#define V_ABS 4
#define SPECIES_ABS 5
#define MESH_DISPLACEMENT1_ABS 6
#define MESH_DISPLACEMENT2_ABS 7
#define MESH_DISPLACEMENT3_ABS 8
#define PRESSURE_ABS 9
#define SHEAR_RATE_ABS 10
#define S11_ABS 11
#define S12_ABS 12
#define S22_ABS 13
#define S11_1_ABS 14
#define S12_1_ABS 15
#define S22_1_ABS 16
#define S13_1_ABS 17
#define S23_1_ABS 18
#define S33_1_ABS 19
#define S11_2_ABS 20
#define S12_2_ABS 21
#define S22_2_ABS 22
#define S13_2_ABS 23
#define S12_3_ABS 24
#define S22_3_ABS 25
#define S13_3_ABS 26
#define S23_3_ABS 27
#define S33_3_ABS 28
#define S11_4_ABS 29
#define S12_4_ABS 30
#define S22_4_ABS 31
#define S13_4_ABS 32
#define S23_4_ABS 33
#define S22_5_ABS 34
#define S13_5_ABS 35
#define S23_5_ABS 36
#define S33_5_ABS 37
#define S11_6_ABS 38
#define S12_6_ABS 39
#define S22_6_ABS 40
#define S13_6_ABS 41
#define S23_6_ABS 42
#define S33_6_ABS 43
#define S13_7_ABS 44
#define S23_7_ABS 45
#define S33_7_ABS 46
#define VELOCITY3_ABS		47
#define TEMPERATURE_ABS		48	
#define MASS_FRACTION_ABS	49	
#define S13_ABS		50	
#define S23_ABS		51	
#define S33_ABS		52	
#define S23_2_ABS	53
#define S11_3_ABS	54
#define S33_4_ABS	55
#define S33_2_ABS	56
#define S11_5_ABS	57
#define S12_5_ABS	58
#define S11_7_ABS	59
#define S12_7_ABS	60
#define S22_7_ABS	61	
	
#define LINEAR_INT      1
#define QUADRATIC_INT   2
#define QUAD_GP_INT     3
#define BIQUADRATIC_INT 4

enum mode {input_file, mat_file, reset};
enum table_types {NO_TABLE, BC_TABLE, VISCOSITY_TABLE, SATURATION_TABLE, SPECIES_TABLE};


#ifndef _RF_IO_H
#define _RF_IO_H

#include "rf_io_const.h"	 /* just in case it has not been done yet.
				  * -> We need this include dor this file */

extern char Input_File[MAX_FNL]; /* input EXODUS II database w/ problem defn */

extern char ExoFile[MAX_FNL];	  /* input EXODUS II database w/ problem defn */

extern char Exo_LB_File[MAX_FNL]; /* EXODUS II load balance info for mesh */

extern char ExoFileOut[MAX_FNL];  /* output EXODUS II database w/ results */

extern char ExoAuxFile[MAX_FNL];  /* auxiliary EXODUS II database for initguess */

extern char DomainMappingFile[MAX_FNL]; /* Domain Mapping file. Maps the materials
				       and names of material boundaries 
				       specified in this file into this file
				       into chemkin domains and chemkin
				       surface and volumetric domains. */

extern int CPU_word_size;

extern int IO_word_size;        /* Precision variables for exodus II files*/

extern char Init_GuessFile[MAX_FNL];
			        /* ASCII file holding initial guess */

extern char Soln_OutFile[MAX_FNL]; 
			        /* ASCII file holding solution vector, */
			        /* same format as Init_GuessFile      */

extern int Debug_Flag;		/* Flag to specify debug info is to be     */
				/* printed out. The value of this flag     */
				/* determines the level of diagnostic info */
				/* which is printed to stdout              */
				/* Debug_Flag == 0 	No output          */
			   	/*	 1	minimun output             */
				/*	 2	medium  output		   */
				/*	 3	maximum output             */
			   	/*	 -1	check jacobian             */
				/*	 -2	check jacobian with scaling*/
#ifdef MATRIX_DUMP
extern int Number_Jac_Dump;     /* Number of jacobians to dump out
				 * If the value is negative, then the one
				 * jacobian, the -n'th jacobian, is dumped
				 * out */
#endif
extern int Iout;		/* Flag to specify level of diagnositc output */
				/* which is to be printed out for the program */

extern int Write_Intermediate_Solutions;
				/* Flag specifies whether to */
			        /* write out solution data at each */
                                /* Newton iteration. */
extern int Write_Initial_Solution; 
				/* Flag to indicate whether to write the
			         * initial solution to the ascii and 
				 * exodus output files */
extern int Num_Var_Init ;	/* Number of variables to overwrite with
				 * global initialization */
extern int Num_Var_External;	/* Number of external variables */
extern int Num_Var_External_pix;/* number of external variables (pixel only)*/
extern int Anneal_Mesh;         /* Flag specifying creation of a special exodus
				 * file with coordinates adjusted to the
				 * deformed coordinates (i.e. new displacements
				 * are set to zero and mesh is deemed stress-free */

extern const char anneal_file[];

extern char front_scratch_directory[MAX_FNL]; /* def'd in main.c */

extern int  Unlimited_Output; /* flag to unlimit multiprocessor output */
                              /*    def'd in main.c    */
#endif


