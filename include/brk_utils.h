#ifndef BRK_UTILS_H
#define BRK_UTILS_H

#include "exo_struct.h"

extern void check_for_brkfile(char* brkfile_name);

extern void write_brk_file(char* brkfile_name, Exo_DB *exo);

extern void call_brk();

extern void fix_output();

#endif
