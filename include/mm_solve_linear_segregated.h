#ifndef _MM_SOLVE_LINEAR_SEGREGATED_H
#define _MM_SOLVE_LINEAR_SEGREGATED_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_SOLVE_LINEAR_SEGREGATED_C
#define EXTERN
#endif

#ifndef _MM_SOLVE_LINEAR_SEGREGATED_C
#define EXTERN extern
#endif

EXTERN int
solve_linear_segregated(struct Aztec_Linear_Solver_System *ams,
			/* ptrs to Aztec linear systems */
			double x[],     /* soln vector on this proc */
			double delta_t, /* time step size */
			double theta,   /* parameter to vary time 
					 * integration from 
					 *   explicit (theta = 1) to 
					 *   implicit (theta = 0) */
			double x_old[], /* soln vector @ previous time */
			double x_older[], /* soln vector @ previous, previous time */
			double xdot[],  /* dxdt predicted for new time */
			double xdot_old[], /* dxdt for previous time */
			double resid_vector[],
			double x_update[], 
			double scale[],  /*Scale factor held for modified newton
					  *resolves */
			int *converged, /* whether the Newton iteration
					 * has converged (out) */
			int *nprint,    /* counter for time step number */
			double *glob_var_vals,   /* global variable values */
			/* details about post proc vars */
			double time_value,
			Exo_DB *exo,
			Dpi *dpi,
			Comm_Ex *cx,
			int nt,
			int *time_step_reform);

EXTERN int 
linear_solve(struct Aztec_Linear_Solver_System *ams,
	     double x[],     /* soln vector on this proc */
	     double resid_vector[],
	     double delta_x[],
	     char stringer[]);


#endif /* _MM_SOLVE_LINEAR_SEGREGATED_H */
