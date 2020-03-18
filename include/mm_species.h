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
 
#ifndef GOMA_MM_SPECIES_H
#define GOMA_MM_SPECIES_H


#include "std.h"

struct Variable_Initialization;

extern int normalize_species_fractions (double [], const int);
extern int check_consistent_fraction_vector (struct Variable_Initialization *,
						   int, int, double[]);

extern double wt_from_Xk (const int, const double *, const double *);
extern double wt_from_Yk (const int, const double *, const double *);
extern double wt_from_Ck (const int, const double *, const double *);
extern void Xk_from_Yk (const int, double *, double *, const double *);
extern void Yk_from_Xk (const int, double *, double *, const double *);
extern void Yk_from_Ck(const int, double *, double *, MATRL_PROP_STRUCT *);
extern void Ck_from_Xk(const int, double *, double *, MATRL_PROP_STRUCT *, const double);
extern void Xk_from_Ck(const int, double *, double *);
extern void Ck_from_Yk(const int, double *, double *, MATRL_PROP_STRUCT *, const double);
extern void Ck_from_Dk(const int, double *, double *, MATRL_PROP_STRUCT *);
extern void Dk_from_Ck(const int, double *, double *, MATRL_PROP_STRUCT *);
extern int convert_species_var (int, MATRL_PROP_STRUCT *,
				      int, double *, double );
extern void deriv1_Ck_to_Yk(double *, MATRL_PROP_STRUCT *, double *, double);
extern void assign_species_var_type(const int, const int, const int);
extern void assign_global_species_var_type(const int, const int);
extern void assign_species_prefix(const int, char *);
#endif
