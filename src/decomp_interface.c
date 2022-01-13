#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "decomp_interface.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "rf_io.h"
#include "mm_as_structs.h"
#include "std.h"

#define MAX_DECOMP_COMMAND 2048

void call_decomp(char *exodus_file) {
  char decomp_command[MAX_DECOMP_COMMAND + 1];

  //  int snerr = snprintf(decomp_command, MAX_DECOMP_COMMAND, "bash -c decomp -p %d %s", Num_Proc,
  //  exodus_file);
  int snerr =
      snprintf(decomp_command, MAX_DECOMP_COMMAND,
               "nem_slice -e -S  -l inertial -c -o initial.exoII.nem -m mesh=2 %s", exodus_file);
  if (Debug_Flag) {
    DPRINTF(stdout, "decomp command: %s\n", decomp_command);
  }
  if (snerr < 0 || snerr >= MAX_DECOMP_COMMAND) {
    GOMA_EH(GOMA_ERROR, "Error creating decomp command snprintf");
  }

  int syserr = system(decomp_command);
  if (syserr != 0) {
    GOMA_EH(GOMA_ERROR, "System call failed for decomp command.");
  }

  if (WEXITSTATUS(syserr) == 127) {
    GOMA_EH(GOMA_ERROR, "System call failed, decomp not found, add SEACAS utils to PATH");
  }
}

void decompose_exodus_files(void) {
  int i;

  if (strcmp(ExoAuxFile, "") != 0) {
    if (Debug_Flag) {
      DPRINTF(stdout, "Decomposing exodus file %s\n", ExoAuxFile);
    }
    call_decomp(ExoAuxFile);
  }

  if (efv->Num_external_field != 0) {
    for (i = 0; i < efv->Num_external_field; i++) {
      if (Debug_Flag) {
        DPRINTF(stdout, "Decomposing exodus file %s\n", efv->file_nm[i]);
      }
      call_decomp(efv->file_nm[i]);
    }
  }
  if (Debug_Flag) {
    DPRINTF(stdout, "Decomposing exodus file %s\n", ExoFile);
  }
  call_decomp(ExoFile);
}
void join_exodus_file(char *filename) { GOMA_EH(GOMA_ERROR, "join exodus file not implemented"); }
