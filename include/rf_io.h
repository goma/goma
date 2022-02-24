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

/*
 * These are declarations for the global I/O varibles needed for
 * FEM problem specification
 *
 *
 * $Id: rf_io.h,v 5.3 2007-12-10 20:59:53 hkmoffa Exp $
 *
 *
 * Note: The associated constants are defined in "rf_io_const.h" and
 *       the definitions are defeind in rf_io_structs.h"
 */

#ifndef GOMA_RF_IO_H
#define GOMA_RF_IO_H

#include "rf_io_const.h" /* just in case it has not been done yet.
				  * -> We need this include dor this file */

extern char Input_File[MAX_FNL]; /* input EXODUS II database w/ problem defn */

extern char ExoFile[MAX_FNL]; /* input EXODUS II database w/ problem defn */

extern char Exo_LB_File[MAX_FNL]; /* EXODUS II load balance info for mesh */

extern char ExoFileOut[MAX_FNL]; /* output EXODUS II database w/ results */

extern char ExoFileOutMono[MAX_FNL]; /* output EXODUSS II database without per proc string */

extern char ExoAuxFile[MAX_FNL]; /* auxiliary EXODUS II database for initguess */

extern int ExoTimePlane; /* initial time plane to use */

extern char Echo_Input_File[MAX_FNL]; /* echo of problem def file  */

extern int Decompose_Flag;

extern int Decompose_Type;

extern char *GomaPetscOptions;
extern int GomaPetscOptionsStrLen;

extern char DomainMappingFile[MAX_FNL]; /* Domain Mapping file. Maps the materials
                                       and names of material boundaries
                                       specified in this file into this file
                                       into chemkin domains and chemkin
                                       surface and volumetric domains. */

extern int CPU_word_size;

extern int IO_word_size; /* Precision variables for exodus II files*/

extern char Init_GuessFile[MAX_FNL];
/* ASCII file holding initial guess */

extern char Soln_OutFile[MAX_FNL];
/* ASCII file holding solution vector, */
/* same format as Init_GuessFile      */

extern int Decompose_Flag; /* Flag to check for built-in brking */

extern int Debug_Flag; /* Flag to specify debug info is to be     */
                       /* printed out. The value of this flag     */
                       /* determines the level of diagnostic info */
                       /* which is printed to stdout              */
                       /* Debug_Flag == 0 	No output          */
                       /*	 1	minimun output             */
                       /*	 2	medium  output		   */
                       /*	 3	maximum output             */
                       /*	 -1	check jacobian             */
                       /*	 -2	check jacobian with scaling*/

extern int New_Parser_Flag; /* New_Parser_Flag = 0	Parse with old parser */
                            /* New_Parser_Flag = 1  Parse with new flex/bison parser */

#ifdef MATRIX_DUMP
extern int Number_Jac_Dump; /* Number of jacobians to dump out
                             * If the value is negative, then the one
                             * jacobian, the -n'th jacobian, is dumped
                             * out */
#endif
extern int Iout; /* Flag to specify level of diagnositc output */
                 /* which is to be printed out for the program */

extern int Write_Intermediate_Solutions;
/* Flag specifies whether to */
/* write out solution data at each */
/* Newton iteration. */
extern int Write_Initial_Solution;
/* Flag to indicate whether to write the
 * initial solution to the ascii and
 * exodus output files */
extern int Num_Var_Init;         /* Number of variables to overwrite with
                                  * global initialization */
extern int Num_Var_Bound;        /* Number of variables to apply bounds */
extern int Num_Var_LS_Init;      /* number of variables to overwirte with
                                  * level set index initialization */
extern int Num_Var_External;     /* Number of external variables */
extern int Num_Var_External_pix; /* number of external variables (pixel only)*/
extern int Anneal_Mesh;          /* Flag specifying creation of a special exodus
                                  * file with coordinates adjusted to the
                                  * deformed coordinates (i.e. new displacements
                                  * are set to zero and mesh is deemed stress-free */

extern double Porous_liq_inventory; /*global variable for finite-insult boundary
                                     * condition
                                     */
extern double **Spec_source_inventory;
extern double *Spec_source_lumped_mass;
extern const char anneal_file[];

extern char front_scratch_directory[MAX_FNL]; /* def'd in main.c */

extern int Unlimited_Output; /* flag to unlimit multiprocessor output */
                             /*    def'd in main.c    */
#endif
