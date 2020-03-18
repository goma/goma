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

/* Utilities for interfacing with brkfix to brk exodus files */

#include <stdio.h>
#include <string.h>

#include "std.h"
#include "exo_struct.h"
#include "brk_utils.h"
#include "brkfix/brk.h"
#include "brkfix/fix.h"
#include "el_elm_info.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "rd_mesh.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_masks.h"
#include "rf_mp.h"

int Brk_Flag;

void
check_for_brkfile(char* brkfile_name) {
  FILE* brkfile;

  /* Check to see if we can read the brkfile */
  brkfile = fopen(brkfile_name, "r");

  if (brkfile == NULL) {
    Brk_Flag = 2;
  } else {
    /*   read_brkfile(brkfile); */
    fclose(brkfile);
  }

  if (Num_Proc > 1 && Brk_Flag == 2) {
    EH(GOMA_ERROR, "Cannot create brkfile on parallel run, run on only 1 processor");
  }

  return;
}

void
write_brk_file(char* brkfile_name, Exo_DB *exo) {
  int i, row, col, eb, eq;
  FILE* brkfile;

  brkfile = fopen(brkfile_name, "w");

  if (brkfile == NULL) {
    EH(GOMA_ERROR, "Cannot open BRK File for writing");
    return;
  }
  
  fprintf(stderr, "Writing brk file to file: %s\n", brkfile_name);

  /* Nbev how many variables are in use */
  /* loop through materials */
  for (eb = 0; eb < upd->Num_Mat; eb++) {
    int var_in_use[MAX_VARIABLE_TYPES];
    int dof_for_var[MAX_VARIABLE_TYPES];
    int nbev;
    int npe = exo->eb_num_nodes_per_elem[eb];
    int elementType;

    /* Find element type for this element */
    elementType = Elem_Type(exo, eb);

    for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
      var_in_use[i] = 0;
      dof_for_var[i] = 0;
    }

    /* loop through elements */
    fprintf(brkfile, "# Material block id\n");
    fprintf(brkfile, "%d\n", exo->eb_id[eb]);
    
    nbev = 0;

    /* Find EQ variables in use */
    for (eq = 0; eq < MAX_VARIABLE_TYPES; eq++) {
      if (pd_glob[eb]->v[0][eq]) {
        var_in_use[nbev] = eq;
        nbev++;
      }
    }

    /* Find max DOF for each variable */
    for (eq = 0; eq < nbev; eq++) {
      int eqid = var_in_use[eq];
      int dof;
      for (i = 0; i < npe; i++) {
        dof = dof_lnode_interp_type(i, elementType, pd_glob[eb]->i[0][eqid], 1);
        if (dof > dof_for_var[eq]) {
          dof_for_var[eq] = dof;
        }
      }
    }
    fprintf(brkfile, "\n");
    fprintf(brkfile, "# Number of basic eqn variables\n");
    fprintf(brkfile, "%d\n", nbev);

    /* Print equation vars and corresponding settings */
    for (eq = 0; eq < nbev; eq++) {
      int eqid = var_in_use[eq];
      if (eqid == MASS_FRACTION) {
        fprintf(brkfile, "%d %d %d %d\n", eqid, 1, pd_glob[eb]->Num_Species, dof_for_var[eq]);
      } else {
        fprintf(brkfile, "%d %d %d %d\n", eqid, 1, 1, dof_for_var[eq]);
      }
    }

    fprintf(brkfile, "\n");
    fprintf(brkfile, "# Interaction Table\n");

    /* find interactions based on variables in use */
    for (row = 0; row < nbev; row++) {
      for (col = 0; col < nbev; col++) {
        fprintf(brkfile, "%d ", Inter_Mask[pg->imtrx][var_in_use[row]][var_in_use[col]]);
      }
      fprintf(brkfile, "\n");
    }

    fprintf(brkfile, "\n");
    fprintf(brkfile, "# Num of nodes per element\n");
    fprintf(brkfile, "%d\n", npe);

    /* Loop through each node and check degree of freedom for each var*/
    for (i = 0; i < npe; i++) {
      for (eq = 0; eq < nbev; eq++) {
        int eqid = var_in_use[eq];
        int dof;

        /* Use edge as 1 for worst case */
        dof = dof_lnode_interp_type(i, elementType, pd_glob[eb]->i[0][eqid], 1);
        if (dof > 0) {
          fprintf(brkfile, "%d ", eqid);
        } 
      }
      fprintf(brkfile, "\n");
    }
    fprintf(brkfile, "\n");
  }

  fclose(brkfile);
  return;
}

void
call_brk(void)
{
  int i;

  if( strcmp( ExoAuxFile, "" ) != 0 ) {
    if (Debug_Flag) {
      DPRINTF(stdout, "Brking exodus file %s\n", ExoAuxFile);
    }
    brk_exo_file(Num_Proc, Brk_File, ExoAuxFile);
  }

  if( efv->Num_external_field != 0 ) {
    for( i=0; i<efv->Num_external_field; i++ ) {
      if (Debug_Flag) {
        DPRINTF(stdout, "Brking exodus file %s\n", efv->file_nm[i]);
      }
      brk_exo_file(Num_Proc, Brk_File, efv->file_nm[i]);
    }
  }
  if (Debug_Flag) {
    DPRINTF(stdout, "Brking exodus file %s\n", ExoFile);
  }
  brk_exo_file(Num_Proc, Brk_File, ExoFile);
}

void
fix_output(void)
{
  if (Debug_Flag) {
    DPRINTF(stdout, "\nFixing exodus file pieces.\n");
  }
  fix_exo_file(Num_Proc, ExoFileOutMono);
}
