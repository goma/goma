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

	AUXILIARY MISC ROUTINES

	BY IAN GATES

	AUGUST 18 1997

	FIRST VERSION FOR GOMA.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "std.h"
#define GOMA_SL_AUXUTIL_C
#include "sl_auxutil.h"

#ifdef aix
#define MV_CSR   mv_csr
#define MV_MSR   mv_msr
#else
#ifdef hpux
#define MV_CSR   MV_CSR
#define MV_MSR   MV_MSR
#else
#ifdef solaris
#define MV_CSR   MV_CSR
#define MV_MSR   MV_MSR
#else
#ifdef sgi
#define MV_CSR   MV_CSR
#define MV_MSR   MV_MSR
#else
#define MV_CSR   MV_CSR
#define MV_MSR   MV_MSR
#endif
#endif
#endif
#endif

/*

	DYNAMIC MEMORY ALLOCATION ROUTINES

	BY IAN GATES

	AUGUST 18 1997

	FIRST VERSION FOR GOMA.

*/
double
**Dmatrix_birth ( int n,
                  int m )
{
	int i, err=0;
	double **a=NULL;
	if (n*m == 0)
		{ return a; }
	if ((a = calloc(n, sizeof(double*))) == NULL) 
		{ err = 1; }
	for (i=0;i<n;i++) {
		if ((a[i] = calloc(m, sizeof(double))) == NULL)
			{ err = 1; }}
	if (err == 1) {
		puts(" DM_birth(): Insufficient memory.");
		exit(0); }
	return a;
} /* END of routine **Dmatrix_birth */
/**************************************************************************/

int
**Imatrix_birth ( int n,
                  int m )
{
	int i, err=0;
	int **a=NULL;
	if (n*m == 0)
		{ return a; }
	if ((a = calloc(n, sizeof(int*))) == NULL)
		{ err = 1; }
	for (i=0;i<n;i++) {
		if ((a[i] = calloc(m, sizeof(int))) == NULL)
			{ err = 1; }}
	if (err == 1) {
		puts(" IM_birth(): Insufficient memory.");
		exit(0); }
	return a;
} /* END of routine **Imatrix_birth */
/**************************************************************************/

double
*Dvector_birth ( int n )
{
	int err=0;
	double *a=NULL;
	if (n == 0)
		{ return a; }
	if ((a = calloc(n, sizeof(double))) == NULL) 
		{ err = 1; }
	if (err == 1) {
		puts(" DV_birth(): Insufficient memory.");
		exit(0); }
	return a;
} /* END of routine **Dvector_birth */
/**************************************************************************/

int
*Ivector_birth ( int n )
{
	int err=0;
	int *a=NULL;
	if (n == 0)
		{ return a; }
	if ((a = calloc(n, sizeof(int))) == NULL)
		{ err = 1; }
	if (err == 1) {
		puts(" IV_birth(): Insufficient memory.");
		exit(0); }
	return a;
}/* END of routine *Ivector_birth */
/**************************************************************************/

void
Dmatrix_death ( double **a,
                int n,
                int m )
{
	int i;
	for (i=0;i<n;i++) 
		{ free(a[i]); }
	free(a);
}/* END of routine Dmatrix_death */
/**************************************************************************/

void
Imatrix_death ( int **a,
                int n,
                int m )
{
	int i;
	for (i=0;i<n;i++)
		{ free(a[i]); }
	free(a);
}/* END of routine Imatrix_death */
/**************************************************************************/


void
Dvector_death ( double *a,
                int n )
{
	free(a);
}/* END of routine Dvector_death */
/**************************************************************************/

void
Ivector_death ( int *a,
                int n )
{
	free(a);
}/* END of routine Ivector_death */
/**************************************************************************/
/*

	VECTOR MANIPULATION ROUTINES

	BY IAN GATES

	SEPT 28 1997

	FIRST VERSION FOR GOMA.

*/
/*
  L2 NORM OF A VECTOR
*/
double 
nnorm(int n, double *v)
{
  int i;
  double f;
  f = 0.0;
  for (i=0;i<n;i++) { f += SQUARE(v[i]); } 
  return sqrt(f);
}
/*
  DOT PRODUCT
*/
double 
dot_product(int n, double *v, double *w)
{
  int i;
  double f;
  f = 0.0;
  for (i=0;i<n;i++) { f += v[i]*w[i]; }	
  return f;
}
/*
  SET VECTOR TO ZERO (aa = 0)
*/
void 
vzero(int n, double *aa)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = 0.0; }
}
/*
  SET VECTOR TO A CONSTANT (aa = a)
*/
void 
vinit(int n, double *aa, double a)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = a; }
}
/*
  COPY ONE VECTOR TO ANOTHER (aa = b*bb)
*/
void 
vcopy(int n, double *aa, double b, double *bb)
{
  int i;
  if (b == 1.0)
     { for (i=0;i<n;i++) { aa[i] = bb[i]; }}
  else
     { for (i=0;i<n;i++) { aa[i] = b*bb[i]; }}
}
/*
  SUM TWO VECTORS (aa = b*bb+c*cc)
*/
void 
v2sum(int n, double *aa, double b, double *bb, double c, double *cc)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = b*bb[i]+c*cc[i]; }
}
/*
  SUM THREE VECTORS (aa = b*bb+c*cc+d*dd)
*/
void 
v3sum(int n, double *aa, double b, double *bb, double c, double *cc, double d, double *dd)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = b*bb[i]+c*cc[i]+d*dd[i]; }
}
/*
  ADD A VECTOR TO ANOTHER VECTOR (aa += b*bb)
*/
void 
v1add(int n, double *aa, double b, double *bb)
{
  int i;
  for (i=0;i<n;i++) { aa[i] += b*bb[i]; }
}
/*
  ADD TWO VECTORS TO A THIRD VECTOR (aa += b*bb+c*cc)
*/
void 
v2add(int n, double *aa, double b, double *bb, double c, double *cc)
{
  int i;
  for (i=0;i<n;i++) { aa[i] += b*bb[i]+c*cc[i]; }
}
/*
  ADD THREE VECTORS TO FOURTH VECTOR (aa += b*bb+c*cc+d*dd)
*/
void 
v3add(int n, double *aa, double b, double *bb, double c, double *cc, double d, double *dd)
{
  int i;
  for (i=0;i<n;i++) { aa[i] += b*bb[i]+c*cc[i]+d*dd[i]; }
}
/*
  CHANGE SIGN OF VECTOR (aa = -aa)
*/
void 
vchange_sign(int n, double *aa)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = -aa[i]; }
}
/*
  PRODUCT OF VECTOR COMPONENTS (aa *= b*bb)
*/
void 
v2product(int n, double *aa, double b, double *bb)
{
  int i;
  for (i=0;i<n;i++) { aa[i] *= b*bb[i]; }
}
/*
  PRODUCT OF TWO VECTORS (aa = bb*cc)
*/
void 
vc_product(int n, double *aa, double *bb, double *cc)
{
  int i;
  for (i=0;i<n;i++) { aa[i] = bb[i]*cc[i]; }
}
/*
  QUOTIENT OF TWO VECTORS (aa = bb/cc)
*/
void 
vc_quotient(int n, double *aa, double *bb, double *cc)
{
  int i;
  for (i=0;i<n;i++) { aa[i] *= bb[i]/cc[i]; }
}
/*
  MAXIMUM ENTRY IN VECTOR (a = max(aa))
*/
double 
vc_max(int n, double *aa)
{
  int i;double a=aa[0];
  for (i=1;i<n;i++) { if (aa[i] > a) { a = aa[i]; }} 
  return a;
}
/*
  MINIMUM ENTRY IN VECTOR (a = min(aa))
*/
double 
vc_min(int n, double *aa)
{
  int i;double a=aa[0];
  for (i=1;i<n;i++) { if (aa[i] < a) { a = aa[i]; }} 
  return a;
}
/*
  SCALAR - VECTOR PRODUCT (aa = a*aa)
*/
void 
vsproduct(int n, double *aa, double a)
{
  int i;
  for (i=0;i<n;i++) { aa[i] *= a; }
}
/*
  MATRIX-VECTOR PRODUCT (CSR)
*/
void 
MV_CSR(int *n, int *raw, int *col, double *a, double *v, double *w)
{
  int i, j, k;
  double *wk;
  wk = Dvector_birth(*n);
  for (i=0;i<*n;i++) {
    k=0;
    for (j=raw[i];j<raw[i+1];j++) {
      wk[k] = v[col[j]];
      k++; }
    w[i] = dot_product(raw[i+1]-raw[i], &a[raw[i]], &wk[0]); }
  Dvector_death(&wk[0], *n);
}
/*
  MATRIX-VECTOR PRODUCT (MSR)
*/
void 
MV_MSR(int *n, int *ija, double *a, double *v, double *w)
{
  int i, j, k;
  double *wk;
  wk = Dvector_birth(*n);
  for (i=0;i<*n;i++)
    { w[i] = a[i]*v[i]; }
  for (i=0;i<*n;i++) {
    k = 0;
    for (j=ija[i];j<ija[i+1];j++) {
      wk[k] = v[ija[j]];
      k++; }
    w[i] += dot_product(k, &a[ija[i]], &wk[0]); }
  Dvector_death(&wk[0], *n);
}
/* END of file sl_auxutil.c  */
/**************************************************************************/

