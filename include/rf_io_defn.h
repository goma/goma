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
 * These are definitions for the global I/O varibles needed for
 * FEM problem specification (i.e., this file allocates storage and
 * can be used to initialize the value of global variables).
 *
 *
 * $Id: rf_io_defn.h,v 5.3 2007-12-10 20:59:53 hkmoffa Exp $
 *
 *
 * Note: The associated constants are defined in "rf_io_const.h" and
 *       the declarations for structures are in rf_io_structs.h"
 *       Declarations for the variables defined in this file are
 *       included in rf_io.h.
 */

#ifndef _RF_IO_DEFN_H
#define _RF_IO_DEFN_H

#include <limits.h>
#include "std.h"
#include "rf_io_const.h"

char	Input_File[MAX_FNL]="\0";  /* input EXODUS II database w/ problem defn */

char	ExoFile[MAX_FNL]="\0";	   /* input EXODUS II database w/ problem defn */

char	Exo_LB_File[MAX_FNL]="\0"; /* EXODUS II load balance info for mesh */

char	ExoFileOut[MAX_FNL]="\0";  /* output EXODUS II database w/ results */

char	ExoFileOutMono[MAX_FNL]="\0";  /* output EXODUS II database without per proc identifier */

char	ExoAuxFile[MAX_FNL]="\0";  /* auxiliary EXODUS II database for initguess */

int     ExoTimePlane = INT_MAX;    /* Time plane # or continuation # of soln to use as an initguess */

char	Echo_Input_File[MAX_FNL]="\0";	/* echo of problem def file  */

char    Brk_File[MAX_FNL]="\0";        /* input file for brk called as subroutine */

char    DomainMappingFile[MAX_FNL]="\0"; /* Domain Mapping file. Maps the materials
				       and names of material boundaries 
				       specified in this file into this file
				       into chemkin domains and chemkin
				       surface and volumetric domains. */

int     CPU_word_size;
int     IO_word_size;		/*Precision variables for exodus II files*/

char	Init_GuessFile[MAX_FNL];/* ASCII file holding initial guess */

char	Soln_OutFile[MAX_FNL];	/* ASCII file holding solution vector, */
				/* same format as Init_GuessFile      */

int	Debug_Flag;		/* Flag to specify debug info is to be     */
				/* printed out. The value of this flag     */
				/* determines the level of diagnostic info */
				/* which is printed to stdout              */
				/* Debug_Flag == 0 	No output          */
			   	/*	 1	minimun output             */
				/*	 2	medium  output		   */
				/*	 3	maximum output             */
			   	/*	 -1	check jacobian             */
				/*	 -2	check jacobian with scaling*/
				
int	New_Parser_Flag;	/* New_Parser_Flag = 0	Parse with old parser */
				/* New_Parser_Flag = 1  Parse with new flex/bison parser */
					
				
#ifdef MATRIX_DUMP
int     Number_Jac_Dump = 0;    /* Number of jacobians to dump out
				 * If the value is negative, then the one
				 * jacobian, the -n'th jacobian, is dumped
				 * out */
#endif
int	Iout;			/* Flag to specify level of diagnostic output */
				/* which is to be printed out for the program */

int Write_Intermediate_Solutions = FALSE; /* Flag specifies whether to */
	      			          /* write out solution data at each */
                                          /* Newton iteration. */
int     Write_Initial_Solution = FALSE; 
				/* Flag to indicate whether to write the
			         * initial solution to the ascii and exodus
				 * output files */
int     Num_Var_Init ;		/* number of variables to overwrite with
				 * global initialization */
int     Num_Var_LS_Init;        /* number of variables to overwirte with
				 * level set index initialization */
int     Num_Var_External;	/* number of total external variables (exoII or pixel)*/
int     Num_Var_External_pix;	/* number of external variables (pixel only)*/
int     Anneal_Mesh;            /* flag specifying creation of a special exodus
				 *  file with coordinates adjusted to the
				 * deformed coordinates (i.e. new displacements
				 * are set to zero and mesh is deemed stress-free */

double Porous_liq_inventory; /*global variable for finite-insult boundary condition*/
double **Spec_source_inventory; /*global variable for cumulative reacted source */
double *Spec_source_lumped_mass; /*global variable for species lumped mass */

const char anneal_file[] = ANNEAL_FILE_NAME;


/*
 * Benner's frontal solver wants to strcat() onto this directory name,
 * so leave enough space for dirname/lu.123456.0, for example.
 *
 * Here it is defined and initialized.
 */

char front_scratch_directory[MAX_FNL] = FRONT_SCRATCH_DIRECTORY;

/******************************************************************************/
/* end of rf_io_defn.h                                                        */
/******************************************************************************/

#endif
