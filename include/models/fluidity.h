#ifndef GOMA_MODELS_FLUIDITY_H
#define GOMA_MODELS_FLUIDITY_H

#ifdef __cplusplus
extern "C" {
#endif
#include "std.h"
#include "mm_as_structs.h"
#include "mm_viscosity.h"

dbl fluidity_viscosity(int fluidity_species, /* integer associated with conc eqn for bond */
                       VISCOSITY_DEPENDENCE_STRUCT *d_mu);

int fluidity_source(int species_no, struct Species_Conservation_Terms *st, dbl *params);

#ifdef __cplusplus
} // extern "C"
#endif
#endif // GOMA_MODELS_FLUIDITY_H