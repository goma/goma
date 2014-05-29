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
 *	sl_rcm.h
 *
 *	RCM DRIVER DEFINES AND ROUTINE PROTOTYPES
 *
 *	3/21/01 PRS
 *
 */

#ifndef _RCM_H
#define _RCM_H
#endif

#define UPPER_CASE		/* Solaris, HP/UX, Linux, IRIX, janus, ... */

#if defined (_AIX)
#undef UPPER_CASE
#define LOWER_CASE
#endif

/*
 * Two cases to consider here:
 * A) Using AztecOO from the current Trilinos distribution:
 *    Here, the Fortran function in az_reorder.f is called az_rcm().
 * Three things can happen to RCM
 *	(1) RCM -> az_rcm
 *	(2) RCM -> az_rcm_
 *	(3) RCM -> az_rcm__
 * This will be taken to be true when TRILINOS is defined in the makefile.
 *
 * B) Using an older version of Aztec or AztecOO:
 *    Here, the Fortran function in az_reorder.f is called rcm().
 * Three things can happen to RCM
 *	(1) RCM -> rcm
 *	(2) RCM -> rcm_
 *	(3) RCM -> rcm__
 * This will be the case if TRILINOS is not defined, or if AZTEC_2_1 is.
 *
 * Say, (2) is most likely backup default unless specified here.
 */

/* Set AZRCM according to makefile defines */
#ifdef AZRCM
#undef AZRCM
#endif

#if ( defined(TRILINOS) && !defined(AZTEC_2_1) )
#define AZRCM
#endif

#ifdef AZRCM
#define RCM az_rcm_
#else
#define RCM rcm_
#endif

#undef  APPEND_0_UNDERSCORE
#define APPEND_1_UNDERSCORE
#undef  APPEND_2_UNDERSCORE

#if defined (_AIX) || defined (__hpux)  
#define APPEND_0_UNDERSCORE
#undef  APPEND_1_UNDERSCORE
#undef  APPEND_2_UNDERSCORE
#endif

#if ( defined(__sun) && defined(__SVR4) ) || defined (SGI) || defined (__PUMAGON__ ) || ( defined (linux) && defined (COMPILER_64BIT) ) 
#undef  APPEND_0_UNDERSCORE
#define APPEND_1_UNDERSCORE
#undef  APPEND_2_UNDERSCORE
#endif

#if ( ( defined (linux) && !defined (COMPILER_64BIT) ) || defined(darwin) ) && !defined(__INTEL_COMPILER) 
#undef  APPEND_0_UNDERSCORE
#undef  APPEND_2_UNDERSCORE
#define APPEND_1_UNDERSCORE
#endif

#ifdef AZRCM  /* Current AztecOO */

#ifdef APPEND_0_UNDERSCORE
#define RCM az_rcm
#endif

#ifdef APPEND_1_UNDERSCORE
#define RCM az_rcm_
#endif

#ifdef APPEND_2_UNDERSCORE
#define RCM az_rcm__
#endif

#else  /* Older AztecOO or Aztec */
#ifdef APPEND_0_UNDERSCORE
#define RCM rcm
#endif

#ifdef APPEND_1_UNDERSCORE
#define RCM rcm_
#endif

#ifdef APPEND_2_UNDERSCORE
#define RCM rcm__
#endif

#endif

