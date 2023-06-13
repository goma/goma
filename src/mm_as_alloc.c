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
 *$Id: mm_as_alloc.c,v 5.19 2010-04-07 22:27:00 prschun Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "el_elm.h"  /* Has shape, element type stuff for bf_init */
#include "el_geom.h" /* Has info I'd like to replicate into the */
#include "rf_allo.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "std.h"
/* Problem_Description structure... */
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "rd_mesh.h"
#include "rf_vars_const.h"

#define GOMA_MM_AS_ALLOC_C

static void init_Viscoelastic_Nonmodal(struct Viscoelastic_Nonmodal *);

static void init_Viscoelastic_Constitutive(struct Viscoelastic_Constitutive *);

static void init_Generalized_Newtonian(GEN_NEWT_STRUCT *);

static void init_Elastic_Constitutive(struct Elastic_Constitutive *);

static void init_Viscoplastic_Constitutive(struct Viscoplastic_Constitutive *);

static int shape_list(Exo_DB *);

// Definitions and initializations of global pointers
/*
 * This is the single location where these are defined.
 */

UPD_STRUCT *upd = NULL;
PROBLEM_DESCRIPTION_STRUCT **pd_glob = NULL, *pd = NULL;
PROBLEM_GRAPH_STRUCT *pg = NULL;
struct Element_Indices **ei = NULL;
struct Element_Indices **eiRelated = {NULL};
struct Element_Stiffness_Pointers *esp = NULL;
struct Element_Quality_Metrics *eqm = NULL;
struct Element_Variable_Pointers *esp_old = NULL, *esp_dot = NULL, *esp_dbl_dot = NULL, *evp = NULL;
struct Action_Flags *af = NULL;
BASIS_FUNCTIONS_STRUCT **bf = NULL;
BASIS_FUNCTIONS_STRUCT **bfd = NULL;
BASIS_FUNCTIONS_STRUCT **bfi = NULL;
BASIS_FUNCTIONS_STRUCT **bfex = NULL;
struct Shell_Block **shell_blocks = NULL;
struct Field_Variables *fv = NULL, *fv_sens = NULL;
struct Diet_Field_Variables *fv_dot_dot = NULL, *fv_dot_dot_old = NULL;
struct Diet_Field_Variables *fv_old = NULL, *fv_dot = NULL;
struct Diet_Field_Variables *fv_dot_old = NULL;
struct Constitutive_Relations **cr_glob = NULL, *cr = NULL;
struct Local_Element_Contributions *lec = NULL;
struct Porous_Media_Variables *pmv = NULL, *pmv_old = NULL;
PMV_ML_STRUCT *pmv_ml = NULL, *pmv_ml_old = NULL;
struct Porous_Media_Variables_Hysteresis *pmv_hyst = NULL;
struct Material_Properties *mp = NULL, **mp_glob = NULL, *mp_old = NULL;
GEN_NEWT_STRUCT *gn = NULL;
GEN_NEWT_STRUCT **gn_glob = NULL;
struct Viscoelastic_Constitutive **ve = NULL, ***ve_glob = NULL;
struct Viscoelastic_Nonmodal *vn = NULL, **vn_glob = NULL;
struct Elastic_Constitutive *elc = NULL, **elc_glob = NULL, *elc_rs = NULL, **elc_rs_glob = NULL;
struct Viscoplastic_Constitutive *evpl = NULL, **evpl_glob = NULL;
struct Transient_Information *tran = NULL;
struct Library_IO *libio = NULL;
struct Eigensolver_Info *eigen = NULL;
struct Continuation_Information *cont = NULL;
struct Loca_Input *loca_in = NULL;
struct AC_Information *augc = NULL;
struct HC_Information *huntc = NULL;
struct External_Field_Variables *efv = NULL;
STABILIZATION_PARAMS_STRUCT *Stab = NULL;

struct Lubrication_Auxiliaries *LubAux = NULL;
struct Lubrication_Auxiliaries *LubAux_old = NULL;

/*
 * This is where these are declared. The pointer is loaded with the
 * address of a double zero that will be used when element degree of freedom
 * pointers need to be assigned to point to something harmless.
 */

dbl *p0;

dbl zero = 0.0;

/*
 * In this memory allocation routine, "pd", etc. have one level of indirection
 * greater than their usual meaning. The main driver "main" calls these
 * routines to make sure that "pd", "ei" point to places with enough
 * space to hold the relevent information.
 */

#include "mm_eh.h"

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int pd_alloc(void)

/******************************************************************
 *
 * pd_alloc():
 *     Assign space for and initialize the uniform problem
 * description structure and the problem description structures.
 ******************************************************************/
{
  int status = 0;
  int mn;
  static char yo[] = "pd_alloc";
  upd = alloc_struct_1(struct Uniform_Problem_Description, 1);
  upd->Species_Var_Type = SPECIES_UNDEFINED_FORM;

  pd_glob = (struct Problem_Description **)alloc_ptr_1(MAX_NUMBER_MATLS);
  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    pd_glob[mn] = alloc_struct_1(struct Problem_Description, 1);
    pd_glob[mn]->Species_Var_Type = SPECIES_UNDEFINED_FORM;
  }

  pg = alloc_struct_1(struct Problem_Graph, 1);

  int imtrx;
  for (imtrx = 0; imtrx < MAX_NUM_MATRICES; imtrx++) {
    pg->time_step_control_disabled[imtrx] = 0;
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Problem_Description @ %p has %ld bytes", yo, (void *)pd_glob,
            (long int)sizeof(struct Problem_Description));
  }
  if (upd == NULL) {
    status = -1;
    GOMA_EH(status, "Hay un problema conseguiendo la memoria para Uniform_Problem_Description");
  }
  if (pd_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Problem_Description");
  }
  if (pg == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Problem_Graph");
  }

  return (status);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int efv_alloc(void) {
  int sz;
  int status = 0;
  static char yo[] = "efv_alloc";

  /*
   * External_Field_Variables___________________________________________________
   */

  sz = sizeof(struct External_Field_Variables);
  efv = (struct External_Field_Variables *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: External_Field_Variables @ %p has %d bytes", yo, (void *)efv, sz);
  }

  if (efv == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for External_Field_Variables");
  }

  memset(efv, 0, sz);

  return (status);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int mp_alloc(void)

/***********************************************************************
 *
 * mp_alloc()
 *
 *  Allocate Material_Properties structures and also allocate
 *  substructures underneath the main structure (there are none as
 *  of yet).
 *
 *  Right now, we allocate a fixed number of material property
 *  structures given by the constant, MAX_NUMBER_MATLS, plus an
 *  additional Material_Property structure, mp_old, which I don't know
 *  the reason for.
 *
 *  Return
 *       0 -> Everything went well
 *       fatal error -> Can't malloc enough memory
 ***********************************************************************/
{
  int mn;
  static char yo[] = "mp_alloc";

  mp_glob = (MATRL_PROP_STRUCT **)alloc_ptr_1(MAX_NUMBER_MATLS);
  if (mp_glob == NULL) {
    GOMA_EH(-1, "Problem getting memory for Material_Properties");
  }

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    mp_glob[mn] = (MATRL_PROP_STRUCT *)alloc_void_struct_1(sizeof(MATRL_PROP_STRUCT), 1);
  }
  mp_old = (MATRL_PROP_STRUCT *)alloc_void_struct_1(sizeof(MATRL_PROP_STRUCT), 1);

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Material_Properties @ %p has %lu bytes", yo, (void *)mp,
            (long unsigned int)sizeof(MATRL_PROP_STRUCT));
  }
  return 0;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int cr_alloc(void)

/******************************************************************
 *
 * cr_alloc:
 *
 *  Initializes the global structure, cr_glob.
 ******************************************************************/
{
  int sz;
  int mn;
  int status;
  static char yo[] = "cr_alloc";

  status = 0;

  /*
   * Constitutive_Relations____________________________________________________
   */

  sz = sizeof(struct Constitutive_Relations *);
  cr_glob = (struct Constitutive_Relations **)array_alloc(1, MAX_NUMBER_MATLS, sz);

  sz = sizeof(struct Constitutive_Relations);

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    cr_glob[mn] = (struct Constitutive_Relations *)array_alloc(1, 1, sz);
    memset(cr_glob[mn], 0, sz);
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Constitutive_Relations @ %p has %d bytes", yo, (void *)cr, sz);
  }

  if (cr_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Constitutive_Relations");
  }

  return (status);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int gn_alloc(void) {
  int sz;
  int mn;
  int status;
  static char yo[] = "gn_alloc";

  status = 0;

  /*
   * Generalized_Newtonian______________________________________________________
   */

  sz = sizeof(struct Generalized_Newtonian *);
  gn_glob = (struct Generalized_Newtonian **)array_alloc(1, MAX_NUMBER_MATLS, sz);

  if (gn_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Generalized_Newtonian");
  }

  sz = sizeof(struct Generalized_Newtonian);

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    gn_glob[mn] = (struct Generalized_Newtonian *)array_alloc(1, 1, sz);

    if (gn_glob[mn] == NULL) {
      status = -1;
      GOMA_EH(status, "Problem getting memory for Generalized_Newtonian material");
    }
    memset(gn_glob[mn], 0, sz);
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Generalized_Newtonian @ %p has %d bytes", yo, (void *)gn, sz);
  }

  return (status);
}

int ve_alloc(void) {
  int sz;
  int mn;
  int mode;
  int status;
  static char yo[] = "ve_alloc";

  status = 0;

  /*
   * Viscoelastic_Constitutive_________________________________________________
   */

  sz = sizeof(struct Viscoelastic_Nonmodal *);
  vn_glob = (struct Viscoelastic_Nonmodal **)smalloc(MAX_NUMBER_MATLS * sz);

  sz = sizeof(struct Viscoelastic_Nonmodal);

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    vn_glob[mn] = (struct Viscoelastic_Nonmodal *)smalloc(sz);
    init_Viscoelastic_Nonmodal(vn_glob[mn]);
  }

  sz = sizeof(struct Viscoelastic_Constitutive *);
  ve = (struct Viscoelastic_Constitutive **)smalloc(MAX_MODES * sz);

  sz = sizeof(struct Viscoelastic_Constitutive *);
  ve_glob = (struct Viscoelastic_Constitutive ***)array_alloc(2, MAX_NUMBER_MATLS, MAX_MODES, sz);

  sz = sizeof(struct Viscoelastic_Constitutive);

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    for (mode = 0; mode < MAX_MODES; mode++) {
      ve_glob[mn][mode] = (struct Viscoelastic_Constitutive *)smalloc(sz);
      init_Viscoelastic_Constitutive(ve_glob[mn][mode]);
    }
  }

  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    for (mode = 0; mode < MAX_MODES; mode++) {
      ve_glob[mn][mode]->gn = alloc_struct_1(GEN_NEWT_STRUCT, 1);
      init_Generalized_Newtonian(ve_glob[mn][mode]->gn);
    }
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Viscoelastic_Constitutive @ %p has %d bytes", yo, (void *)ve_glob, sz);
  }

  if (ve_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Viscoelastic_Constitutive");
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int elc_alloc(void) {
  int sz;
  int mn;
  int status;
  static char yo[] = "elc_alloc";

  status = 0;

  /*
   * Elastic_Constitutive______________________________________________________
   */

  sz = sizeof(struct Elastic_Constitutive *);
  elc_glob = (struct Elastic_Constitutive **)smalloc(MAX_NUMBER_MATLS * sz);

  sz = sizeof(struct Elastic_Constitutive);
  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    elc_glob[mn] = (struct Elastic_Constitutive *)smalloc(sz);
    memset(elc_glob[mn], 0, sz);
    init_Elastic_Constitutive(elc_glob[mn]);
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Elastic Constitutive @ %p has %d bytes", yo, (void *)gn_glob, sz);
  }

  if (elc_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Elastic_Constitutive");
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int evp_alloc(void) {
  int sz;
  int mn;
  int status;
  static char yo[] = "evp_alloc";

  status = 0;

  /*
   * Viscoplastic_Constitutive_________________________________________________
   */

  sz = sizeof(struct Viscoplastic_Constitutive *);
  evpl_glob = (struct Viscoplastic_Constitutive **)smalloc(MAX_NUMBER_MATLS * sz);

  sz = sizeof(struct Viscoplastic_Constitutive);
  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    evpl_glob[mn] = (struct Viscoplastic_Constitutive *)smalloc(sz);
    init_Viscoplastic_Constitutive(evpl_glob[mn]);
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Viscoplastic Constitutive @ %p has %d bytes", yo, (void *)gn_glob, sz);
  }

  if (evpl_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Viscoplastic_Constitutive");
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int evp_tensor_alloc(Exo_DB *exo) {
  int status;
  int ielem_type, ip_total, ielem0, mn;
  int i, ip, a, b, c, k, w;
  status = 0;

  /*
   * Viscoplastic global tensors_______________________________________________
   */

  for (mn = 0; mn < upd->Num_Mat; mn++) {
    /*Don't bother with allocation if we are not solving an EVP model */
    if (evpl_glob[mn]->ConstitutiveEquation == EVP_HYPER) {
      pd_glob[mn]->Num_Dim = Num_Dim; /* from "el_geom.h" */
      ielem0 = exo->eb_ptr[mn];
      ielem_type = Elem_Type(exo, ielem0);     /* element type */
      ip_total = elem_info(NQUAD, ielem_type); /* number of */
      /* quadrature points */

      if (evpl_glob[mn] == NULL) {
        status = -1;
        GOMA_EH(status, "Problem getting memory for Viscoplastic_Constitutive");
      }

      evpl_glob[mn]->F_vp_glob =
          (dbl ****)array_alloc(4, exo->num_elems, ip_total, DIM, DIM, sizeof(dbl));
      evpl_glob[mn]->F_vp_old_glob =
          (dbl ****)array_alloc(4, exo->num_elems, ip_total, DIM, DIM, sizeof(dbl));
      evpl_glob[mn]->TT_glob =
          (dbl ****)array_alloc(4, exo->num_elems, ip_total, DIM, DIM, sizeof(dbl));
      evpl_glob[mn]->TT_old_glob =
          (dbl ****)array_alloc(4, exo->num_elems, ip_total, DIM, DIM, sizeof(dbl));
      evpl_glob[mn]->dTT_dx_glob =
          (dbl ******)array_alloc(6, exo->num_elems, ip_total, DIM, DIM, DIM, MDE, sizeof(dbl));
      evpl_glob[mn]->dTT_dx_old_glob =
          (dbl ******)array_alloc(6, exo->num_elems, ip_total, DIM, DIM, DIM, MDE, sizeof(dbl));
      evpl_glob[mn]->dTT_dc_glob = (dbl ******)array_alloc(6, exo->num_elems, ip_total, DIM, DIM,
                                                           MAX_CONC, MDE, sizeof(dbl));
      evpl_glob[mn]->dTT_dc_old_glob = (dbl ******)array_alloc(6, exo->num_elems, ip_total, DIM,
                                                               DIM, MAX_CONC, MDE, sizeof(dbl));

      for (i = 0; i < exo->num_elems; i++) {
        for (ip = 0; ip < ip_total; ip++) {
          for (a = 0; a < DIM; a++) {
            for (b = 0; b < DIM; b++) {
              evpl_glob[mn]->F_vp_old_glob[i][ip][a][b] = 0.;
              evpl_glob[mn]->F_vp_glob[i][ip][a][b] = 0.;
              evpl_glob[mn]->TT_glob[i][ip][a][b] = 0.;
              evpl_glob[mn]->TT_old_glob[i][ip][a][b] = 0.;

              for (k = 0; k < MDE; k++) {
                for (c = 0; c < DIM; c++) {
                  evpl_glob[mn]->dTT_dx_glob[i][ip][a][b][c][k] = 0.;
                  evpl_glob[mn]->dTT_dx_old_glob[i][ip][a][b][c][k] = 0.;
                }
                for (w = 0; w < MAX_CONC; w++) {
                  evpl_glob[mn]->dTT_dc_glob[i][ip][a][b][w][k] = 0.;
                  evpl_glob[mn]->dTT_dc_old_glob[i][ip][a][b][w][k] = 0.;
                }
              }
            }
          }
        }
      }
    }
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int elc_rs_alloc(void) {
  int sz;
  int mn;
  int status;
  static char yo[] = "elc_alloc";

  status = 0;

  /*
   * Elastic_Constitutive____________________________________________________
   */

  sz = sizeof(struct Elastic_Constitutive *);
  elc_rs_glob = (struct Elastic_Constitutive **)smalloc(MAX_NUMBER_MATLS * sz);

  sz = sizeof(struct Elastic_Constitutive);
  for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
    elc_rs_glob[mn] = (struct Elastic_Constitutive *)smalloc(sz);
    init_Elastic_Constitutive(elc_rs_glob[mn]);
  }

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Elastic Constitutive RS @ %p has %d bytes", yo, (void *)gn_glob, sz);
  }

  if (elc_rs_glob == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Elastic_Constitutive");
  }

  return (status);
}

int tran_alloc(void) {
  int sz;
  int status;
  static char yo[] = "tran_alloc";

  status = 0;

  /*
   * Transient_Information______________________________________________________
   */

  sz = sizeof(struct Transient_Information);
  tran = (struct Transient_Information *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Transient_Information @ %p has %d bytes", yo, (void *)tran, sz);
  }

  if (tran == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Transient_Information");
  }

  /*
   * Element_Quality_Metrics___________________________________________________
   */
  sz = sizeof(struct Element_Quality_Metrics);
  eqm = alloc_struct_1(struct Element_Quality_Metrics, 1);

  if (Debug_Flag) {
    P0PRINTF("%s: Element_Quality_Metrics @ %p has %d bytes", yo, (void *)eqm, sz);
  }

  if (eqm == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Element_Quality_Metrics");
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int libio_alloc(void) {
  int sz;
  int status;
  static char yo[] = "libio_alloc";

  status = 0;

  /*
   * Library_IO
   */

  sz = sizeof(struct Library_IO);
  libio = (struct Library_IO *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Library_IO @ %p has %d bytes", yo, (void *)libio, sz);
  }

  if (libio == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Library_IO");
  }

  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int eigen_alloc(void) {
  int sz;
  int status;
  static char yo[] = "eigen_alloc";

  status = 0;

  /*
   * Eigensolver_Info
   */

  sz = sizeof(struct Eigensolver_Info);
  eigen = (struct Eigensolver_Info *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Eigensolver_Info @ %p has %d bytes", yo, (void *)eigen, sz);
  }

  if (eigen == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Eigensolver_Info");
  }

  return (status);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

int cont_alloc(void) {
  int sz;
  int status;
  static char yo[] = "cont_alloc";

  status = 0;

  /*
   * Continuation_Information
   */

  sz = sizeof(struct Continuation_Information);
  cont = (struct Continuation_Information *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Continuation_Information @ %p has %d bytes", yo, (void *)cont, sz);
  }

  if (cont == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Continuation_Information");
  }

  return (status);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

int loca_alloc(void) {
  int sz;
  int status;
  static char yo[] = "loca_alloc";

  status = 0;

  /*
   * Loca_Input
   */

  sz = sizeof(struct Loca_Input);
  loca_in = (struct Loca_Input *)array_alloc(1, 1, sz);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Loca_Input @ %p has %d bytes", yo, (void *)cont, sz);
  }

  if (loca_in == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Loca_Input");
  }

  return (status);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
static void Element_Indices_alloc(struct Element_Indices *ei_ptr) {

  ei_ptr->Num_Lvdesc_Per_Var_Type = alloc_int_1(MAX_VARIABLE_TYPES + MAX_CONC, 0);
  ei_ptr->Lvdesc_First_Var_Type = alloc_int_1(MAX_VARIABLE_TYPES + MAX_CONC, -1);
  ei_ptr->Lvdesc_to_ledof = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_to_lvdof = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_to_Var_Type = alloc_int_1(MAX_LOCAL_VAR_DESC, -1);
  ei_ptr->Lvdesc_to_MatID = alloc_int_1(MAX_LOCAL_VAR_DESC, -1);
  ei_ptr->Lvdesc_Numdof = alloc_int_1(MAX_LOCAL_VAR_DESC, 0);
  ei_ptr->Lvdesc_to_Gnn = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_to_Gun = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_to_MFSubvar = alloc_int_1(MAX_LOCAL_VAR_DESC, 0);
  ei_ptr->Lvdesc_to_Lnn = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_Lnn_to_Offset = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_Lnn_to_lvdof = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->lvdof_to_lvdesc = alloc_int_2(MAX_VARIABLE_TYPES + MAX_CONC, MDE, -1);
  ei_ptr->lvdof_to_lvdesc_dof = alloc_int_2(MAX_VARIABLE_TYPES + MAX_CONC, MDE, -1);
  ei_ptr->Lvdesc_Lnn_Numdof = alloc_int_2(MAX_LOCAL_VAR_DESC, MDE, -1);
  ei_ptr->Lvdesc_vd_ptr = (VARIABLE_DESCRIPTION_STRUCT **)alloc_ptr_1(MAX_LOCAL_VAR_DESC);
  /*
   * At this point we don't know the value of Num_Var_Info_Records
   * so we can't malloc the following:
   *       ei_ptr->VDindex_to_Lvdesc = alloc_int_1(Num_Var_Info_Records, -1);
   */

  ei_ptr->owningElementForColVar = alloc_int_1(MAX_VARIABLE_TYPES + MAX_CONC, 0);
  /*
   * For each equation/variable, tell how many dofs must be accounted for
   * in this element and their names(nodes)...
   *
   */
  ei_ptr->dof = alloc_int_1(MAX_VARIABLE_TYPES, 0);
  ei_ptr->dof_ext = alloc_int_1(MAX_EXTERNAL_FIELD, 0);

  ei_ptr->Baby_Dolphin = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);

  /*
   * In more detail, tell the list of nodes associated with those dofs
   * for each variable.
   */
  ei_ptr->dof_list = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);

  /*
   * Finally, pointers into the global unknown arrays for each dof in this
   * element...
   */
  ei_ptr->gnn_list = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);

  /*
   * Finally, pointers into the global unknown arrays for each dof in this
   * element...
   */
  ei_ptr->gun_list = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);
  ei_ptr->ln_to_dof = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);
  ei_ptr->ln_to_first_dof = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);
  ei_ptr->lvdof_to_row_lvdof = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);
  ei_ptr->lvdof_to_ledof = alloc_int_2(MAX_VARIABLE_TYPES, MDE, -1);

  ei_ptr->owned_ledof = alloc_int_1(MDE * MAX_PROB_VAR, 1);
  ei_ptr->ieqn_ledof = alloc_int_1(MDE * MAX_PROB_VAR, -1);
  ei_ptr->matID_ledof = alloc_int_1(MDE * MAX_PROB_VAR, -1);
  ei_ptr->active_interp_ledof = alloc_int_1(MDE * MAX_PROB_VAR, -1);

  ei_ptr->MFsubvar_Offset = alloc_int_2(MAX_CONC, MDE, 0);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int assembly_alloc(Exo_DB *exo)

/********************************************************************
 *
 * assembly_alloc:
 *
 *
 ********************************************************************/
{
  int vi, sz;
  int type; /* index for unique basis functions */
  int status = 0;
  int interp; /* index for interpolation order */
  int si;     /* index for element shape */
  int mn;
  int ipore;
  int num_species_eqn; /* active number of species eqn */
  int imtrx;

  /*
   * The problem description has already been set up. But we need to access
   * it here...
   */

  static char yo[] = "assembly_alloc";

  /*
   * These are critical for the element dof pointers...
   */
  p0 = &zero;

  if (Debug_Flag) {
    DPRINTF(stderr, "%s...\n", yo);
  }

  /*
   * Element_Indices___________________________________________________________
   */

  ei = malloc(sizeof(struct Element_Indices *) * upd->Total_Num_Matrices);
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    ei[pg->imtrx] = alloc_struct_1(struct Element_Indices, 1);
    Element_Indices_alloc(ei[pg->imtrx]);
  }
  pg->imtrx = 0;

  eiRelated = (struct Element_Indices **)alloc_ptr_1(MAX_ELEMENT_INDICES_RELATED);
  for (mn = 0; mn < MAX_ELEMENT_INDICES_RELATED; mn++) {
    eiRelated[mn] = alloc_struct_1(struct Element_Indices, 1);
    Element_Indices_alloc(eiRelated[mn]);
  }

  /*
   * Element_Variable_Pointers________________________________________________
   */

  sz = sizeof(struct Element_Variable_Pointers);
  esp_old = (struct Element_Variable_Pointers *)array_alloc(1, 1, sz);
  esp_dot = (struct Element_Variable_Pointers *)array_alloc(1, 1, sz);
  esp_dbl_dot = (struct Element_Variable_Pointers *)array_alloc(1, 1, sz);
  evp = (struct Element_Variable_Pointers *)array_alloc(1, 1, sz);

  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Element_Variable_Pointers @ %p has %d bytes", yo, (void *)esp_old, sz);
  }

  if ((esp_old == NULL) || (esp_dot == NULL)) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Element_Variable_Pointers");
  }

  /*
   * Element_Stiffness_Pointers________________________________________________
   */
  sz = sizeof(struct Element_Stiffness_Pointers);
  esp = alloc_struct_1(struct Element_Stiffness_Pointers, 1);

  if (Debug_Flag) {
    P0PRINTF("%s: Element_Stiffness_Pointers @ %p has %d bytes", yo, (void *)esp, sz);
  }

  /*
   *  num_species_eqn is equal to the maximum number of species equations
   *  in any one domain
   */
  num_species_eqn = upd->Max_Num_Species_Eqn;

  sz = MDE;

  /* MMH I have included some coupling between particle velocity and
   * other things.  I have skipped over a bunch (like energy) b/c there
   * are no models, yet, that include both energy/temperature and
   * particle velocity
   */

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    /* ENERGY */
    if (Num_Var_In_Type[imtrx][TEMPERATURE]) {
      esp->T = (dbl **)alloc_ptr_1(MDE);
    }

    /* MOMENTUM  */
    if (Num_Var_In_Type[imtrx][VELOCITY1]) {
      esp->v = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    if (Num_Var_In_Type[imtrx][USTAR]) {
      esp->v_star = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* MMH
     * Some of these are just leftovers from copying the velocity stuff.  I left
     * them here in case their particle equivalents are incorporated later.
     */
    /* PARTICLE MOMENTUM */
    if (Num_Var_In_Type[imtrx][PVELOCITY1]) {
      esp->pv = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* MESH_DISPLACEMENT  */
    if (Num_Var_In_Type[imtrx][MESH_DISPLACEMENT1]) {
      esp->d = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* SOLID_DISPLACEMENT  */
    if (Num_Var_In_Type[imtrx][SOLID_DISPLACEMENT1]) {
      esp->d_rs = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* SPECIES CONTINUITY */
    if (Num_Var_In_Type[imtrx][MASS_FRACTION]) {
      esp->c = (dbl ***)alloc_ptr_2(num_species_eqn, MDE);
    }

    /*CONTINUITY */
    if (Num_Var_In_Type[imtrx][PRESSURE]) {
      esp->P = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][PSTAR]) {
      esp->P_star = (dbl **)alloc_ptr_1(MDE);
    }

    /* POLYMER STRESS for all modes */
    if (Num_Var_In_Type[imtrx][POLYMER_STRESS11]) {
      esp->S = (dbl *****)array_alloc(4, MAX_MODES, VIM, VIM, MDE, sizeof(dbl *));
      (void)memset(esp->S[0][0][0], 0, MAX_MODES * VIM * VIM * MDE * sizeof(dbl *));
    }

    /* VELOCITY_GRADIENT */
    if (Num_Var_In_Type[imtrx][VELOCITY_GRADIENT11]) {
      esp->G = (dbl ****)array_alloc(3, VIM, VIM, MDE, sizeof(dbl *));
      (void)memset(esp->G[0][0], 0, VIM * VIM * MDE * sizeof(dbl *));
    }

    /* POTENTIAL */
    if (Num_Var_In_Type[imtrx][VOLTAGE]) {
      esp->V = (dbl **)alloc_ptr_1(MDE);
    }

    /* SURF_CHARGE */
    if (Num_Var_In_Type[imtrx][SURF_CHARGE]) {
      esp->qs = (dbl **)alloc_ptr_1(MDE);
    }

    /* FILL */
    if (Num_Var_In_Type[imtrx][FILL]) {
      esp->F = (dbl **)alloc_ptr_1(MDE);
    }

    /* SHEAR_RATE INVARIANT */
    if (Num_Var_In_Type[imtrx][SHEAR_RATE]) {
      esp->SH = (dbl **)alloc_ptr_1(MDE);
    }

    /* ENORM, |E| */
    if (Num_Var_In_Type[imtrx][ENORM]) {
      esp->Enorm = (dbl **)alloc_ptr_1(MDE);
    }

    /* LEVEL SET CURVATURE */
    if (Num_Var_In_Type[imtrx][CURVATURE]) {
      esp->H = (dbl **)alloc_ptr_1(MDE);
    }

    /* LEVEL SET NORMAL VECTOR */
    if (Num_Var_In_Type[imtrx][NORMAL1]) {
      esp->n = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* POROUS MEDIA VARS */
    if (Num_Var_In_Type[imtrx][POR_LIQ_PRES]) {
      esp->p_liq = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][POR_GAS_PRES]) {
      esp->p_gas = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][POR_POROSITY]) {
      esp->porosity = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][POR_TEMP]) {
      esp->T = (dbl **)alloc_ptr_1(MDE);
    }

    /* VORTICITY PRINCIPLE FLOW DIRECTION
     * This is always 3D */
    if (Num_Var_In_Type[imtrx][VORT_DIR1]) {
      esp->vd = (dbl ***)alloc_ptr_2(DIM, MDE);
    }

    /* Lagrange Multipliers */
    if (Num_Var_In_Type[imtrx][LAGR_MULT1]) {
      esp->lm = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* Structural Shell equations */
    if (Num_Var_In_Type[imtrx][SHELL_CURVATURE]) {
      esp->sh_K = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_CURVATURE2]) {
      esp->sh_K2 = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_TENSION]) {
      esp->sh_tens = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_X]) {
      esp->sh_x = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_Y]) {
      esp->sh_y = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_USER]) {
      esp->sh_u = (dbl **)alloc_ptr_1(MDE);
    }

    /* Shell orientation angles */
    if (Num_Var_In_Type[imtrx][SHELL_ANGLE1]) {
      esp->sh_ang = (dbl ***)alloc_ptr_2(DIM - 1, MDE);
    }

    /*Sundry pieces for Surface Rheological Const. Eqn. */
    if (Num_Var_In_Type[imtrx][R_SHELL_SURF_DIV_V]) {
      esp->div_s_v = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][R_SHELL_SURF_CURV]) {
      esp->curv = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][R_N_DOT_CURL_V]) {
      esp->n_dot_curl_s_v = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][R_GRAD_S_V_DOT_N1]) {
      esp->grad_v_dot_n = (dbl ***)alloc_ptr_2(DIM, MDE);
    }

    if (Num_Var_In_Type[imtrx][R_SHELL_DIFF_FLUX]) {
      esp->sh_J = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][R_SHELL_DIFF_CURVATURE]) {
      esp->sh_Kd = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][R_SHELL_NORMAL1]) {
      esp->n = (dbl ***)alloc_ptr_2(DIM, MDE);
    }

    /* EIGENVALUE FOR VORTICITY PRINCIPLE FLOW DIRECTION
     * This is always 3D */
    if (Num_Var_In_Type[imtrx][VORT_LAMBDA]) {
      esp->vlambda = (dbl **)alloc_ptr_1(MDE);
    }

    /* BOND Evolution
     * associated with structure formation during flow */
    if (Num_Var_In_Type[imtrx][BOND_EVOLUTION]) {
      esp->nn = (dbl **)alloc_ptr_1(MDE);
    }

    /* Extension Velocity */
    if (Num_Var_In_Type[imtrx][EXT_VELOCITY]) {
      esp->ext_v = (dbl **)alloc_ptr_1(MDE);
    }

    /* Electric Field */
    if (Num_Var_In_Type[imtrx][EFIELD1]) {
      esp->E_field = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    /* Phase function */
    if (Num_Var_In_Type[imtrx][PHASE1]) {
      esp->pF = (dbl ***)alloc_ptr_2(pfd->num_phase_funcs, MDE);
    }
    /* Acoustic pressure */
    if (Num_Var_In_Type[imtrx][ACOUS_PREAL]) {
      esp->apr = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][ACOUS_PIMAG]) {
      esp->api = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][ACOUS_REYN_STRESS]) {
      esp->ars = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][EM_CONT_REAL]) {
      esp->epr = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][EM_CONT_IMAG]) {
      esp->epi = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_BDYVELO]) {
      esp->sh_bv = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_LUBP]) {
      esp->sh_p = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][LUBP]) {
      esp->lubp = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][LUBP_2]) {
      esp->lubp_2 = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_FILMP]) {
      esp->sh_fp = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_FILMH]) {
      esp->sh_fh = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_PARTC]) {
      esp->sh_pc = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_SAT_CLOSED]) {
      esp->sh_sat_closed = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_PRESS_OPEN]) {
      esp->sh_p_open = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_PRESS_OPEN_2]) {
      esp->sh_p_open_2 = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_TEMPERATURE]) {
      esp->sh_t = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_DELTAH]) {
      esp->sh_dh = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_LUB_CURV]) {
      esp->sh_l_curv = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_LUB_CURV_2]) {
      esp->sh_l_curv_2 = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_SAT_GASN]) {
      esp->sh_sat_gasn = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][POR_SINK_MASS]) {
      esp->sink_mass = (dbl **)alloc_ptr_1(MDE);
    }
    /* Poynting Vector  */
    if (Num_Var_In_Type[imtrx][LIGHT_INTP] || Num_Var_In_Type[imtrx][LIGHT_INTM] ||
        Num_Var_In_Type[imtrx][LIGHT_INTD]) {
      esp->poynt = (dbl ***)alloc_ptr_2(VIM, MDE);
    }
    /* EM_wave components  */
    if (Num_Var_In_Type[imtrx][EM_E1_REAL] || Num_Var_In_Type[imtrx][EM_E2_REAL] ||
        Num_Var_In_Type[imtrx][EM_E3_REAL]) {
      esp->em_er = (dbl ***)alloc_ptr_2(VIM, MDE);
    }
    if (Num_Var_In_Type[imtrx][EM_E1_IMAG] || Num_Var_In_Type[imtrx][EM_E2_IMAG] ||
        Num_Var_In_Type[imtrx][EM_E3_IMAG]) {
      esp->em_ei = (dbl ***)alloc_ptr_2(VIM, MDE);
    }
    if (Num_Var_In_Type[imtrx][EM_H1_REAL] || Num_Var_In_Type[imtrx][EM_H2_REAL] ||
        Num_Var_In_Type[imtrx][EM_H3_REAL]) {
      esp->em_hr = (dbl ***)alloc_ptr_2(VIM, MDE);
    }
    if (Num_Var_In_Type[imtrx][EM_H1_IMAG] || Num_Var_In_Type[imtrx][EM_H2_IMAG] ||
        Num_Var_In_Type[imtrx][EM_H3_IMAG]) {
      esp->em_hi = (dbl ***)alloc_ptr_2(VIM, MDE);
    }

    if (Num_Var_In_Type[imtrx][TFMP_PRES]) {
      esp->tfmp_pres = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][TFMP_SAT]) {
      esp->tfmp_sat = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][RESTIME]) {
      esp->restime = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][MOMENT0]) {
      esp->moment = (dbl ***)alloc_ptr_2(MAX_MOMENTS, MDE);
    }

    if (Num_Var_In_Type[imtrx][DENSITY_EQN]) {
      esp->rho = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_SHEAR_TOP]) {
      esp->sh_shear_top = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_SHEAR_BOT]) {
      esp->sh_shear_bot = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_CROSS_SHEAR]) {
      esp->sh_cross_shear = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][MAX_STRAIN]) {
      esp->max_strain = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][CUR_STRAIN]) {
      esp->cur_strain = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][SHELL_SAT_1]) {
      esp->sh_sat_1 = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_SAT_2]) {
      esp->sh_sat_2 = (dbl **)alloc_ptr_1(MDE);
    }
    if (Num_Var_In_Type[imtrx][SHELL_SAT_3]) {
      esp->sh_sat_3 = (dbl **)alloc_ptr_1(MDE);
    }

    if (Num_Var_In_Type[imtrx][EDDY_NU]) {
      esp->eddy_nu = (dbl **)alloc_ptr_1(MDE);
    }

  } /* End of loop over matrices */

  /*
   * Action_Flags______________________________________________________________
   */
  af = alloc_struct_1(struct Action_Flags, 1);
  if (Debug_Flag) {
    P0PRINTF("%s: Action_Flags @ %p has %lu bytes", yo, (void *)af,
             (long unsigned int)sizeof(struct Action_Flags));
  }

  /*
   * Basis functions: allocate space only for a list of unique basis functions
   * required to assemble terms for this problem. Then, depending on the
   * particular weighing and interpolation to be used for each variable,
   * point to the appropriate member of the unique set.
   *
   * Here, the "bfd" point to real stuff.
   *
   * Then, the "bf" just point to the right place (see bf_init)...
   */

  /*
   * Hardwired for now...later, load these up based on the EXODUS II database
   */

  /* This was the old way:
    if ( Num_Dim == 2 )
      {
        Num_Shapes       = 1;
        Unique_Shapes[0] = QUADRILATERAL;
        **      Unique_Shapes[1] = TRIANGLE; **
      }
    else if ( Num_Dim == 3 )
      {
        Num_Shapes       = 1;
        Unique_Shapes[0] = HEXAHEDRON;
        **      Unique_Shapes[1] = TETRAHEDRON; **
        **      Unique_Shapes[2] = PRISM; **
      }
    else if ( Num_Dim == 1 )
      {
        Num_Shapes       = 1;
        Unique_Shapes[0] = LINE_SEGMENT;
      }
    else
      {
        GOMA_EH(GOMA_ERROR, "Cannot classify shapes...");
      }
  */

  /* Now, here is the new way: */
  Num_Shapes = shape_list(exo);

  /*
   * To determine how many types of basis functions we'll need to represent,
   * multiply the number of basic element shapes by the number of different
   * interpolations that have been specified...
   */
  Num_Basis_Functions = Num_Shapes * Num_Interpolations;

  bfd = (struct Basis_Functions **)alloc_ptr_1(Num_Basis_Functions);
  bfi = (struct Basis_Functions **)alloc_ptr_1(MAX_INTERP_TYPES);
  bf = (struct Basis_Functions **)alloc_ptr_1(MAX_VARIABLE_TYPES);
  bfex = (struct Basis_Functions **)alloc_ptr_1(MAX_EXTERNAL_FIELD);

  sz = sizeof(struct Basis_Functions);

  /*
   * Now allocate and do some initialization of the fundamental basis functions
   */

  type = 0;
  for (si = 0; si < Num_Shapes; si++) {
    for (interp = 0; interp < Num_Interpolations; interp++) {
      bfd[type] = alloc_struct_1(struct Basis_Functions, 1);
      bfd[type]->Var_Type_MatID = alloc_int_1(upd->Num_Mat, -1);
      bfd[type]->interpolation = Unique_Interpolations[interp];
      bfd[type]->element_shape = Unique_Shapes[si];
      bfd[type]->Max_Dofs_Interpolation = getdofs(Unique_Shapes[si], Unique_Interpolations[interp]);
      /*
       * Find a representative variable type in each material
       * that employs the current interpolation.
       * Store it in the structure.
       *
       * This is an overly broad interpretation. Specifically,
       * the intersection of an element shape and a material
       * type is never tried.
       */
      for (mn = 0; mn < upd->Num_Mat; mn++) {
        for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
          for (vi = 0; ((bfd[type]->Var_Type_MatID[mn] == -1) && (vi < MAX_VARIABLE_TYPES)); vi++) {
            if (pd_glob[mn]->i[imtrx][vi] == bfd[type]->interpolation) {
              bfd[type]->Var_Type_MatID[mn] = vi;
            }
          }
        }
      }
      type++;
    }
  }

  /* set pointers for basis function interpolation types */
  for (interp = 0; interp < Num_Interpolations; interp++) {
    type = Unique_Interpolations[interp];
    bfi[type] = bfd[interp];
  }
  /* so now bfi[pd_glob[mn]->i[eqn]] points to a basis function with
     the correct interpolation */

  /*
   * Field_Variables___________________________________________________________
   */

  fv = alloc_struct_1(struct Field_Variables, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Field_Variables @ %p has %ld bytes", yo, (void *)fv,
            (long int)sizeof(struct Field_Variables));
  }
  if (fv == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Field_Variables");
  }
  memset(fv, 0, sizeof(struct Field_Variables));

  fv_sens = alloc_struct_1(struct Field_Variables, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Field_Variables Sensitivity @ %p has %lu bytes", yo, (void *)fv_sens,
            sizeof(struct Field_Variables));
  }
  if (fv_sens == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Field_Variables (sens)");
  }
  memset(fv_sens, 0, sizeof(struct Field_Variables));

  sz = sizeof(struct Diet_Field_Variables);
  fv_old = alloc_struct_1(struct Diet_Field_Variables, 1);
  fv_dot = alloc_struct_1(struct Diet_Field_Variables, 1);
  fv_dot_dot = alloc_struct_1(struct Diet_Field_Variables, 1);
  fv_dot_old = alloc_struct_1(struct Diet_Field_Variables, 1);
  fv_dot_dot_old = alloc_struct_1(struct Diet_Field_Variables, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Diet_Field_Variables @ %p has %d bytes", yo, (void *)fv_old, sz);
  }
  if (fv_old == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Diet_Field_Variables");
  }

  /*
   * Porous_Media_Variables__(if needed)_________________________________________
   */
  for (mn = 0; mn < upd->Num_Mat && pmv == NULL; mn++) {
    if (mp_glob[mn]->PorousMediaType == POROUS_SATURATED ||
        mp_glob[mn]->PorousMediaType == POROUS_UNSATURATED ||
        mp_glob[mn]->PorousMediaType == POROUS_TWO_PHASE ||
        mp_glob[mn]->PorousMediaType == POROUS_SHELL_UNSATURATED) {
      pmv = alloc_struct_1(struct Porous_Media_Variables, 1);
      pmv_old = alloc_struct_1(struct Porous_Media_Variables, 1);
      break;
    }
  }
  for (mn = 0; mn < upd->Num_Mat && pmv_ml == NULL; mn++) {
    if (mp_glob[mn]->PorousMediaType == POROUS_SATURATED ||
        mp_glob[mn]->PorousMediaType == POROUS_UNSATURATED ||
        mp_glob[mn]->PorousMediaType == POROUS_TWO_PHASE ||
        mp_glob[mn]->PorousMediaType == POROUS_SHELL_UNSATURATED) {
      if ((mp_glob[mn]->Porous_Mass_Lump) ||
          (mp_glob[mn]->PorousMediaType == POROUS_SHELL_UNSATURATED)) {
        pmv_ml = alloc_struct_1(PMV_ML_STRUCT, 1);
      }
      break;
    }
  }
  for (mn = 0; mn < upd->Num_Mat && pmv_hyst == NULL; mn++) {
    if (mp_glob[mn]->PorousMediaType == POROUS_SHELL_UNSATURATED) {
      for (ipore = 0; ipore < pd_glob[mn]->Num_Porous_Shell_Eqn; ipore++) {
        if ((mp_glob[mn]->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST) ||
            (mp_glob[mn]->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST_EXT)) {
          pmv_hyst = alloc_struct_1(struct Porous_Media_Variables_Hysteresis, 1);
          break;
        }
      }
      if (pmv_hyst != NULL) {
        for (ipore = 0; ipore < pd_glob[mn]->Num_Porous_Shell_Eqn; ipore++) {
          pmv_hyst->curve_type[ipore] = alloc_int_1(exo->num_nodes, -1);
          pmv_hyst->curve_type_old[ipore] = alloc_int_1(exo->num_nodes, -1);
          pmv_hyst->curve_switch[ipore] = alloc_int_1(exo->num_nodes, -1);
          pmv_hyst->num_switch[ipore] = alloc_int_1(exo->num_nodes, -1);

          pmv_hyst->sat_switch[ipore] = alloc_dbl_1(exo->num_nodes, 0.0);
          pmv_hyst->cap_pres_switch[ipore] = alloc_dbl_1(exo->num_nodes, 0.0);
          pmv_hyst->sat_min_imbibe[ipore] = alloc_dbl_1(exo->num_nodes, 0.0);
          pmv_hyst->sat_max_drain[ipore] = alloc_dbl_1(exo->num_nodes, 0.0);
        }
      }
      break;
    }
  }

  /*
   * Stabilization_Variables__(if needed)___________________________________
   */
  for (mn = 0; mn < upd->Num_Mat && Stab == NULL; mn++) {
    mp = mp_glob[mn];
    if (mp->Porous_wt_funcModel == SUPG || vn_glob[mn]->wt_funcModel == SUPG ||
        mp->Spwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG || mp->Ewt_funcModel == SUPG) {
      Stab = alloc_struct_1(STABILIZATION_PARAMS_STRUCT, 1);
      break;
    }
  }

  /*
   * Local_Element_Contributions___________________________________________
   */
  lec = alloc_struct_1(struct Local_Element_Contributions, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Local_Element_Contributions @ %p has %ld bytes", yo, (void *)lec,
            (long int)sizeof(struct Local_Element_Contributions));
  }
  if (lec == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Local_Element_Contributions");
  }

  /*
   * Lubrication Auxiliaries_________________________________________________
   */
  LubAux = alloc_struct_1(struct Lubrication_Auxiliaries, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Lubrication_Auxiliaries @ %p has %ld bytes", yo, (void *)LubAux,
            (long int)sizeof(struct Lubrication_Auxiliaries));
  }
  if (LubAux == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Lubrication_Auxiliaries");
  }
  memset(LubAux, 0, sizeof(struct Lubrication_Auxiliaries));

  LubAux_old = alloc_struct_1(struct Lubrication_Auxiliaries, 1);
  if (Debug_Flag) {
    DPRINTF(stdout, "%s: Lubrication_Auxiliaries Old @ %p has %ld bytes", yo, (void *)LubAux_old,
            (long int)sizeof(struct Lubrication_Auxiliaries));
  }
  if (LubAux_old == NULL) {
    status = -1;
    GOMA_EH(status, "Problem getting memory for Lubrication_Auxiliaries Old");
  }
  memset(LubAux_old, 0, sizeof(struct Lubrication_Auxiliaries));

  return (status);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int bf_init(Exo_DB *exo)

/***********************************************************************
 * bf_init:
 *
 * Initialize some key pieces of each of the unique basis functions "bfd"
 * and set up the bf[variable] to point to the right "bfd"...
 *
 * For now, lets assume that bf[eqn] and bf[var] will be identical
 * references to the same basis function.
 *
 ***********************************************************************/
{
  int ebi, m, t = 0, v, status = 0;
  int imtrx;
  int el, shape, type;

  /*
   * For now, assume variable interpolations
   * and equation weightings are the same.
   *
   * We also assume that the interpolation for each variable is the
   * same for all materials in which it is active.
   *
   * Just make each var/eqn's bf point to the appropriate fundamental basis
   * function...
   */

  for (ebi = 0; ebi < Proc_Num_Elem_Blk; ebi++) {
    m = Matilda[ebi];
    if (m >= 0 && exo->eb_num_elems[ebi] > 0) {
      /* These are needed to check for matching element shapes */
      el = exo->eb_ptr[ebi];
      type = Elem_Type(exo, el);
      shape = type2shape(type);

      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
          /*
           * Is this variable active in the problem?
           */
          if (pd_glob[m]->v[imtrx][v]) {
            /*
             * If so, then check to see which prototype basis function of
             * bfd[] is its match...
             * NEW: Check both interpolation AND element shape!
             */
            for (t = 0; t < Num_Basis_Functions; t++) {
              if (pd_glob[m]->i[imtrx][v] == bfd[t]->interpolation &&
                  shape == bfd[t]->element_shape) {
                if (bf[v] != 0) {
                  if (bf[v] != bfd[t]) {
                    GOMA_WH(-1, "bf[] struct isn't general enough to handle"
                                " variables with multiple interpolations in a single problem\n"
                                " We will carefully reset bf for each element");
                  }
                }
                bf[v] = bfd[t];
              }
            }
            if (bf[v] == NULL) {
              GOMA_EH(-1, "Could not find a match for variable.");
            }
          }
        }
      }
    }
  }

  return (status);
}

/*
 * Initialize some key pieces of each of the unique basis functions "bfd"
 * and set up the bf[variable] to point to the right "bfd"...
 *
 * For now, lets assume that bf[eqn] and bf[var] will be identical references
 * to the same basis function.
 *
 * Material specific matchup done.
 */
int bf_mp_init(struct Problem_Description *pd) {
  int ifound;
  int t, v;
  int status;
  int ishape;

  status = 0;

  t = 0;

  /* This is needed to check for matching element shapes */
  ishape = ei[pg->imtrx]->ielem_shape;

  /*
   * For now, assume variable interpolations
   * and equation weightings are the same.
   *
   * Just make each var/eqn's bf point to the appropriate fundamental basis
   * function...
   */
  for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
    /*
     * Is this variable active in the material?
     *  If it isn't, then don't overwrite the entry, because we seem
     *  to need an idea for the interpolation even if its wrong in order
     *  for internal-boundary problems to work correctly.
     */
    int imtrx;
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (pd->v[imtrx][v]) {
        /*
         * If so, then check to see which prototype basis function of
         * bfd[] is its match...
         * Check both interpolation AND element shape!
         */
        bf[v] = NULL;
        for (t = 0; t < Num_Basis_Functions; t++) {
          if ((pd->i[imtrx][v] == bfd[t]->interpolation) && (ishape == bfd[t]->element_shape))

            if ((pd->i[imtrx][v] == bfd[t]->interpolation) && (ishape == bfd[t]->element_shape)) {
              bf[v] = bfd[t];
            }
        }
        if (bf[v] == NULL) {
          GOMA_EH(-1, "Could not find a match for variable.");
        }
      }
    }
  }
  /*
   * Repeat the same proceedure for external fixed field interpolations
   */

  /*
   * Are there external fields defined in the problem?
   */
  if (efv->ev) {
    for (v = 0; v < efv->Num_external_field; v++) {
      /*
       * If so, then check to see which prototype basis function of
       * bfd[] is its match...
       */
      ifound = 0;
      for (t = 0; t < Num_Basis_Functions; t++) {
        if (efv->i[v] == bfd[t]->interpolation) {
          bfex[v] = bfd[t];
          ifound = 1;
        }
      }
      if (!ifound && efv->i[v] != I_TABLE) {
        GOMA_EH(GOMA_ERROR, "Could not find a match for EXTERNAL variable. Match some active field "
                            "interpolation with those of external variables");
      }
      if (bfex[v] == NULL && efv->i[v] != I_TABLE) {
        GOMA_EH(-1, "Could not find a match for EXTERNAL variable.");
      }
    }
  }

  return (status);
}

/*
 * Initialize the basis function representation to the NULL vector
 */
void bf_reset(void) {
  int v;
  for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
    bf[v] = NULL;
  }
  if (efv->ev) {
    for (v = 0; v < efv->Num_external_field; v++) {
      bfex[v] = NULL;
    }
  }
}

/*
 * Initialize a freshly allocated structure to some safe defaults.
 *
 * Created: 2000/01/27 14:21 MST pasacki@sandia.gov
 *
 * Revised:
 */
static void init_Viscoelastic_Nonmodal(struct Viscoelastic_Nonmodal *v) {
  v->ConstitutiveEquation = 0;
  v->wt_func = (double)0;
  v->wt_funcModel = 0;
  v->eps = (double)0;
  v->evssModel = 0;
  v->modes = 0;
  v->dg_J_model = 0;
  v->dg_J_model_wt = NULL;
  v->len_dg_J_model_wt = 0;
  v->shiftModel = 0;
  v->shift = NULL;
  v->len_shift = 0;
  return;
}

/*
 * Initialize a freshly allocated structure to some safe defaults.
 *
 * Created: 2000/01/27 14:34 MST pasacki@sandia.gov
 *
 * Revised:
 */

static void init_Viscoelastic_Constitutive(struct Viscoelastic_Constitutive *v) {
  v->gn = NULL;
  v->time_const = (double)0;
  v->time_constModel = 0;
  v->alpha = (double)0;
  v->alphaModel = 0;
  v->xi = (double)0;
  v->xiModel = 0;
  v->eps = (double)0;
  v->epsModel = 0;
  v->pos_ls.alpha = 0.0;
  v->pos_ls.eps = 0.0;
  v->pos_ls.xi = 0.0;
  v->pos_ls.time_const = 0.0;
  return;
}

/*
 * Initialize a freshly allocated structure to some safe defaults.
 *
 * Created: 2000/01/27 14:38 MST pasacki@sandia.gov
 *
 * Revised:
 */

static void init_Generalized_Newtonian(struct Generalized_Newtonian *g) {
  g->ConstitutiveEquation = 0;
  g->mu0 = (double)0;
  g->mu0Model = 0;
  g->nexp = (double)0;
  g->nexpModel = 0;
  g->muinf = (double)0;
  g->muinfModel = 0;
  g->lam = (double)0;
  g->lamModel = 0;
  g->aexp = (double)0;
  g->aexpModel = 0;
  g->atexp = (double)0;
  g->atexpModel = 0;
  g->wlfc2 = (double)0;
  g->wlfc2Model = 0;
  g->tau_y = (double)0;
  g->tau_yModel = 0;
  g->fexp = (double)0;
  g->fexpModel = 0;
  g->maxpack = (double)0;
  g->maxpackModel = 0;
  g->sus_species_no = 0;
  g->gelpoint = (double)0;
  g->gelpointModel = 0;
  g->cureaexp = (double)0;
  g->cureaexpModel = 0;
  g->curebexp = (double)0;
  g->curebexpModel = 0;
  g->tgel0 = (double)0;
  g->tgel0Model = 0;
  g->cure_species_no = 0;
  g->DilViscModel = 0;
  g->DilVisc0 = 0.0;

  return;
}

/*
 * Initialize a freshly allocated structure to some safe defaults.
 *
 * Created: 2000/01/27 14:49 MST pasacki@sandia.gov
 *
 * Revised:
 */

static void init_Elastic_Constitutive(struct Elastic_Constitutive *e) {
  int i;

  e->ConstitutiveEquation = 0;
  e->lame_mu = (double)0;
  e->lame_mu_model = 0;
  e->len_u_mu = 0;
  e->u_mu = NULL;
  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    e->d_lame_mu[i] = (double)0;
  }
  e->lame_mu_tableid = 0;

  e->lame_lambda = (double)0;
  e->lame_lambda_model = 0;
  e->len_u_lambda = 0;
  e->u_lambda = NULL;
  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    e->d_lame_lambda[i] = (double)0;
  }

  e->lame_TempShift = 0;
  e->lameTempShiftModel = 0;
  e->len_u_lame_TempShift = 0;
  e->u_lame_TempShift = NULL;
  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    e->d_lame_TempShift[i] = 0;
  }
  e->lame_TempShift_tableid = 0;

  e->bend_stiffness = 0;
  e->bend_stiffness_model = 0;
  e->len_u_bend_stiffness = 0;
  e->u_bend_stiffness = NULL;

  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    e->d_bend_stiffness[i] = 0;
  }

  e->poisson = (double)0;
  e->Strss_fr_sol_vol_frac = (double)0;

  for (i = 0; i < DIM; i++) {
    e->v_mesh_sfs[i] = (double)0;
  }

  e->v_mesh_sfs_model = 0;
  e->len_u_v_mesh_sfs = 0;
  e->u_v_mesh_sfs = NULL;

  e->thermal_expansion = (double)0;
  e->thermal_expansion_model = 0;
  e->len_u_thermal_expansion = 0;
  e->u_thermal_expansion = NULL;
  e->solid_reference_temp = (double)0;
  e->solid_reference_temp_model = 0;

  return;
}

static void init_Viscoplastic_Constitutive(struct Viscoplastic_Constitutive *v) {
  int i;

  v->ConstitutiveEquation = 0;
  v->update_flag = 0;

  v->plastic_mu = (double)0;
  v->plastic_mu_model = 0;
  v->len_u_plastic_mu = 0;
  v->u_plastic_mu = NULL;
  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    v->d_plastic_mu[i] = (double)0;
  }
  v->yield = (double)0;
  v->yield_model = 0;
  v->len_u_yield = 0;
  v->u_yield = NULL;
  for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
    v->d_yield[i] = (double)0;
  }
  v->F_vp_glob = NULL;
  v->F_vp_old_glob = NULL;
  v->TT_glob = NULL;
  v->TT_old_glob = NULL;
  v->dTT_dx_glob = NULL;
  v->dTT_dx_old_glob = NULL;
  v->dTT_dc_glob = NULL;
  v->dTT_dc_old_glob = NULL;
  return;
}

static int shape_list(Exo_DB *exo)
/*
 * Determines number and type of element shapes in use.
 * Sets Num_Shapes and populates Unique_Shapes array.
 */
{
  int e, e1, e2, eshape, etype, there, N;

  /* Elements to search */
  e1 = exo->eb_ptr[0];
  e2 = exo->eb_ptr[exo->num_elem_blocks];

  /* Initialize shape list with first element */
  etype = Elem_Type(exo, e1);
  eshape = type2shape(etype);
  Unique_Shapes[0] = eshape;
  N = 1;

  /* Loop over all other elements, check shapes */
  for (e = e1 + 1; e < e2; e++) {
    etype = Elem_Type(exo, e);
    eshape = type2shape(etype);
    there = in_list(eshape, 0, N, Unique_Shapes);

    /* Add to list if not there */
    if (there == -1) {
      Unique_Shapes[N] = eshape;
      N++;
    }
  }

  if (N > 1)
    DPRINTF(stdout, "\nFound %d different element shapes.\n", N);
  return N;
}

/* end of file mm_as_alloc.c */
