#include "Simmath.h"

#ifdef __cplusplus 
#include <iostream>
using std::printf;
#endif


#define PROBLEM_DIMENSION 2


/* adapted from itkLBFGSOptimizerTest.cxx */
main() {

    smHandle optimizer;
    double results[PROBLEM_DIMENSION];
    int i;
    double params[10];
    void costFunc( int, double*, double*, double*, void* );

    optimizer = smMallocOptimizer( PROBLEM_DIMENSION);

    params[0] = 0;
    smSetOptimizerParameters( optimizer, TRACE, params );

    params[0] = 100;
    smSetOptimizerParameters( optimizer, MAX_FUNCTION_EVALUATIONS, params );

    params[0] =  .0001;
    smSetOptimizerParameters( optimizer, GRADIENT_CONVERGENCE_TOLERANCE, params );

    params[0] =  1.0;
     smSetOptimizerParameters( optimizer, DEFAULT_STEP_LENGTH, params);

    params[0] =  0.9;
    smSetOptimizerParameters( optimizer, LINE_SEARCH_ACCURACY, params);

    smSetCostFunction( optimizer, costFunc );

    // We start not so far from  | 2 -2 |

    results[0] =  100;
    results[1] = -100;

    smDumpOptimizerState(optimizer);

    smRunOptimizer( optimizer, results );
    for( i=0; i<PROBLEM_DIMENSION; i++ ) {
       printf(" results[%d] = %f \n",i,results[i]); 
    }

    smFreeOptimizer(optimizer);
}
   
void costFunc( int n, double *position, double *f, double *g, void* user_data ) {
  double x, y;

   x = position[0]; 
   y = position[1];  

   f[0] = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y; 
   g[0] = 3*x + 2*y -2;
   g[1] = 2*x + 6*y +8; 

   return;
}
