#ifndef GOMA_OMEGA_H_INTERFACE_H
#define GOMA_OMEGA_H_INTERFACE_H

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

typedef struct Exodus_Database Exo_DB;
typedef struct Distributed_Processing_Information Dpi;

void adapt_mesh_omega_h(struct Aztec_Linear_Solver_System **ams,
                        Exo_DB *exo,
                        Dpi *dpi,
                        double **x,
                        double **x_old,
                        double **x_older,
                        double **xdot,
                        double **xdot_old,
                        double **x_oldest,
                        double **resid_vector,
                        double **x_update,
                        double **scale,
                        int step);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif
#endif

// vim: expandtab sw=2 ts=8
