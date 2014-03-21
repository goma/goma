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
 
#ifndef _MM_SPECIES_H
#define _MM_SPECIES_H


extern int normalize_species_fractions PROTO((double [], const int));
extern int check_consistent_fraction_vector PROTO((struct Variable_Initialization *,
						   int, int, double[]));

extern double wt_from_Xk PROTO((const int, const double *, const double *));
extern double wt_from_Yk PROTO((const int, const double *, const double *));
extern double wt_from_Ck PROTO((const int, const double *, const double *));
extern void Xk_from_Yk PROTO((const int, double *, double *, const double *));
extern void Yk_from_Xk PROTO((const int, double *, double *, const double *));
extern void Yk_from_Ck(const int, double *, double *, MATRL_PROP_STRUCT *);
extern void Ck_from_Xk(const int, double *, double *, MATRL_PROP_STRUCT *, const double);
extern void Xk_from_Ck(const int, double *, double *);
extern void Ck_from_Yk(const int, double *, double *, MATRL_PROP_STRUCT *, const double);
extern int convert_species_var PROTO((int, MATRL_PROP_STRUCT *,
				      int, double *, double ));
extern void deriv1_Ck_to_Yk(double *, MATRL_PROP_STRUCT *, double *, double);
extern void assign_species_var_type(const int, const int, const int);
extern void assign_global_species_var_type(const int, const int);
extern void assign_species_prefix(const int, char *);
#endif
