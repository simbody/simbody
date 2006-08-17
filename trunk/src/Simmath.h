#ifndef _SimTK_SIMMATH_H_
#define _SimTK_SIMMATH_H_

typedef void* smHandle; 
typedef int  smStatus;

enum { 	OPTIMIZER, ROOT_FINDER };
enum { LBFGS, POWELL, GRADIENT_DESENT };
enum { SUCCESS, MALLOC_FAILED, UNKNOWN_OPTIMIZER, UNKNOWN_PARAMETER, INVALID_VALUE };
enum { INITIAL_VALUES, TRACE, MAX_FUNCTION_EVALUATIONS, DEFAULT_STEP_LENGTH, LINE_SEARCH_ACCURACY, GRADIENT_CONVERGENCE_TOLERANCE };

#ifdef __cplusplus
namespace SimTK {
class smObject {
   public:
    void *data;
};
} //  namespace SimTK

extern "C" {
#endif

extern smHandle smMallocOptimizer(int, smStatus*);
extern smStatus smDumpOptimizerState( smHandle);
extern smStatus smSetOptimizerParameters( smHandle, unsigned int, double*);
extern smStatus smGetOptimizerParameters( smHandle, unsigned int, double*);
extern smStatus smSetCostFunction( smHandle, void (*costFunction)(double*,double*,double*) );
extern smStatus smRunOptimizer( smHandle, double * );
extern void     smFreeOptimizer( smHandle );

#ifdef __cplusplus
}  /* extern "C" */

#endif


#endif //_SimTK_SIMMATH_H_
