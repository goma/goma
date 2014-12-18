#ifndef INCLUDE_SL_STRATIMIKOS_INTERFACE_H_
#define INCLUDE_SL_STRATIMIKOS_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

int stratimikos_solve(struct Aztec_Linear_Solver_System *ams, double *x_,
    double *b_, int *iterations, char *stratimikos_file);

#ifdef __cplusplus
} // end of extern "C"
#endif

#endif /* INCLUDE_SL_STRATIMIKOS_INTERFACE_H_ */
