#include "Simmath.h"

#ifdef __cplusplus 
#include <iostream>
using std::printf;
#endif


#define PROBLEM_DIMENSION 2


/* adapted from itkLBFGSOptimizerTest.cxx */
main() {

    smHandle optimizer;
    smStatus status;
    double results[PROBLEM_DIMENSION];
    int i;
    double params[10];
    void costFunc( double*, double*, double* );

    optimizer = smMallocOptimizer( PROBLEM_DIMENSION, &status);

    if( status != SUCCESS ){
       printf("malloc failed status = %d\n", status); 
       exit(0);
    }

    params[0] = 0;
    status = smSetOptimizerParameters( optimizer, TRACE, params );
    if( status != SUCCESS ) printf("set TRACE failed status = %d\n", status);

    params[0] = 100;
    status = smSetOptimizerParameters( optimizer, MAX_FUNCTION_EVALUATIONS, params );
    if( status != SUCCESS ) printf("set MAX_FUNCTION_EVALUATIONS failed status = %d\n", status);

    params[0] =  .0001;
    status = smSetOptimizerParameters( optimizer, GRADIENT_CONVERGENCE_TOLERANCE, params );
    if( status != SUCCESS ) printf("set GRADIENT_CONVERGENCE_TOLERANCE failed status = %d\n", status);

    params[0] =  1.0;
    status = smSetOptimizerParameters( optimizer, DEFAULT_STEP_LENGTH, params);
    if( status != SUCCESS ) printf("set GRADIENT_CONVERGENCE_TOLERANCE failed status = %d\n", status);

    params[0] =  0.9;
    status = smSetOptimizerParameters( optimizer, LINE_SEARCH_ACCURACY, params);
    if( status != SUCCESS ) printf("set LINE_SEARCH_ACCURACY failed status = %d\n", status);

    status = smSetCostFunction( optimizer, costFunc );
    if( status != SUCCESS ) printf("set cost function failed status = %d\n", status);

    // We start not so far from  | 2 -2 |

    results[0] =  100;
    results[1] = -100;

    status = smDumpOptimizerState(optimizer);
    if( status != SUCCESS ) printf("Dump Optimizer state failed status = %d\n", status);

    status = smRunOptimizer( optimizer, results );
    if( status != SUCCESS )  {
       printf("Run Optimizer failed status = %d\n", status);
    }else {
       for( i=0; i<PROBLEM_DIMENSION; i++ ) {
          printf(" results[%d] = %f \n",i,results[i]); 
       }
    }

    smFreeOptimizer(optimizer);

}

   
void costFunc( double *position, double *f, double *g ) {

  int i;

   double x, y;

   x = position[0]; 
   y = position[1];  

   f[0] = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y; 
   g[0] = 3*x + 2*y -2;
   g[1] = 2*x + 6*y +8; 


   return;
}
