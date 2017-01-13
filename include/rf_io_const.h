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
 *$Id: rf_io_const.h,v 5.4 2010-04-07 17:48:19 sarober Exp $
 */

#ifndef _RF_IO_CONST_H
#define _RF_IO_CONST_H

#ifndef MAX_NNV
#define	MAX_NNV		100	/* maximum number of nodal variables */
#endif

#ifndef MAX_NEV
#define	MAX_NEV		25	/* maximum number of element variables */
#endif

#ifndef MAX_NHV
#define	MAX_NHV		 1	/* maximum number of history variables */
#endif

#ifndef MAX_NGV
#define	MAX_NGV		300	/* maximum number of global variables */
#endif

#ifndef MAX_FNL
#define	MAX_FNL		128	/* maximum filename len (EXODUS II, Chemkin)*/
#endif

#ifndef NEXO
#define	NEXO		 2	/* number EXODUS II databases available. */
				/* Convention here: */
				/*      [0] == input EXODUS II db, read-only */
				/*      [1] == output EXODUS II db, writable */
#endif

#ifndef MAX_CHAR_IN_INPUT
#define MAX_CHAR_IN_INPUT	256 /* How many characters in input line? */
#endif

#ifndef MAX_INFO
#define MAX_INFO	101	/* maximum number of "info" records */
#endif

#ifndef MAX_QA
#define MAX_QA		20	/* maximum number of QA records */
#endif

#ifndef MAX_VAR_NAME_LNGTH
#define MAX_VAR_NAME_LNGTH 20 /* maximum length of variable names         */
                              /* for Exodus II db output                  */
                              /* HKM -> Changed it from 10 to             */
                              /*         20 to conform to exodus standard */
                              /*        Note: Chemkin needs at least 16   */
                              /*              Plus 3 for a prefix standard */
#endif

#ifndef MAX_NUMBER_PARAMS
#define MAX_NUMBER_PARAMS      100 /* maximum number of parameters allowed on */
				   /* input cards  */
#endif

#ifndef MAX_VAR_DESC_LNGTH
#define MAX_VAR_DESC_LNGTH	80 /* maximum length of variable descriptions */
				   /* for Exodus II db output */
#endif

/* Definitions of constants for command line options */
#ifndef MAX_COMMAND_LINE_LENGTH
#define MAX_COMMAND_LINE_LENGTH	(1024) /* maximum number of characters in */
				       /* a command line sent to system() */
#endif

#define INPUT_FILE               1 /* command for reading from alternate input file */
#define NO_DISPLAY               2 /* command for eliminating output to screen */
#define APREPRO                  3 /* command to system call aprepro */
#define INEXOII_FILE             4 /* redirect in.exoII file */
#define OUTEXOII_FILE            5 /* redirect out.exoII file */
#define CONTIN_FILE              6 /* redirect contin.dat file */
#define SOLN_FILE                7 /* redirect soln.dat file */
#define DEBUG_OPTION             8 /* change debug flag */
#define FASTQ                    9 /* command to system call fastq */
#define FASTQ_APREPRO           10 /* command to system call aprepro on fastq file and call fastq */
#define EX1EX2V2                11 /* command to system call ex1ex2v2 */
#define BLOT                    12 /* command to system call blot after run is finished */
#define STDOUT_OUT              13 /* command to redirect stdout to file */
#define STDERR_OUT              13 /* command to redirect stderr to file */
#define RELAXATION              14 /* change overrelaxation parameter */
#define HELP_USAGE		15 /* print helpful usage of options, args */
#define NEWTON_NUM              16 /* change debug flag */
#define PRINT_CODE_VERSION	17 /* do exactly that */
#define WORKING_DIRECTORY	18 /* working directory to run in */
#define PARSER_OPTION		19 /* command to use new flex/bison parser */
#define NOECHO              20 /* Disable echoing of input/material files */
#define TIME_START              21 /* Initial simulation time */
#define TIME_END                22 /* Maximum simulation time */

#define CONT_BEG_PVALUE        101 /* BEGIN VALUE */
#define CONT_END_PVALUE        102 /* END VALUE */
#define CONT_PATHSTEP          103 /* PATH STEP */
#define CONT_PATH_STEPS        104 /* NUMBER OF PATH STEPS */
#define CONT_MIN_PVALUE        105 /* MINIMUM PATH STEP SIZE */
#define CONT_MAX_PVALUE        106 /* MAXIMUM PATH STEP SIZE */
#define CONT_METHOD            111 /* METHOD; 0th, 1st, ... */
#define CONT_TYPE              112 /* TYPE; BC(1), MAT(2) */
#define CONT_BCID              121 /* BCID */
#define CONT_DFID              122 /* DATA FLOAT ID */
#define CONT_MPID              131 /* MPID */
#define CONT_MTID              132 /* MAT PROPOERTY TAG ID */
#define CONT_BC_LIST           900 /* BC LIST AND STOP */
#define WRITE_INTERMEDIATE     901 /* Turn Write_Intermediate_Solution on */
#define EXOII_TIME_PLANE       902 /* Specify read_exoII_file time plane (step number) */

#ifndef ANNEAL_FILE_NAME
#define ANNEAL_FILE_NAME		"anneal.exoII"
#endif

/* #ifndef P_tmpdir
   #define p_tmpdir			"/var/tmp"
   #endif  */

#ifndef p_tmpdir
#define p_tmpdir			"/tmp"
#endif

#ifndef FRONT_SCRATCH_DIRECTORY
#define FRONT_SCRATCH_DIRECTORY		p_tmpdir
#endif

#define DP_PROC_PRINT_LIMIT      4 /*  maximum number of processors for   */
                                   /*  which certain kinds of information */
                                   /*  files will be generated during     */
                                   /*  parallel processing runs           */

#endif

