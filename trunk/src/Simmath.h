#ifndef _SimTK_SIMMATH_H_
#define _SimTK_SIMMATH_H_

typedef void* smHandle; 
typedef int  smStatus;

enum { LBFGS, POWELL, GRADIENT_DESENT };
enum { SUCCESS, MALLOC_FAILED, UNKNOWN_OPTIMIZER, UNKNOWN_PARAMETER, INVALID_VALUE };
enum { INITIAL_VALUES, TRACE, MAX_FUNCTION_EVALUATIONS, DEFAULT_STEP_LENGTH, LINE_SEARCH_ACCURACY, GRADIENT_CONVERGENCE_TOLERANCE };

#ifdef __cplusplus
extern "C" {
#endif

extern smHandle smMallocOptimizer(int, int, smStatus*);
extern smStatus smGetOptimizerResults( smHandle, double*);
extern smStatus smDumpOptimizerState( smHandle);
extern smStatus smSetOptimizerParameters( smHandle, unsigned int, double*);
extern smStatus smSetCostFunction( smHandle, void (*costFunction)(double*,double*,double*) );
extern smStatus smStartOptimizer( smHandle );
extern smStatus smFreeOptimizer( smHandle );


#ifdef __cplusplus
}  /* extern "C" */
#endif


#endif //_SimTK_SIMMATH_H_
