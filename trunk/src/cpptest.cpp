#include "Simmath.h"
#include "optimizer.h"

#include <iostream>
using std::cout;
using std::endl;


#define PROBLEM_DIMENSION 2


/* adapted from itkLBFGSOptimizerTest.cxx */
main() {

    double params[10];
    double results[10];
    int status,i;
    void costFunc( double*, double*, double* );

    cout << "cpptest " << endl;
    SimTK::smOptimizer *opt = new SimTK::smOptimizer( PROBLEM_DIMENSION );

    params[0] = 0;
    status = opt->setOptimizerParameters( TRACE, params );
    if( status != SUCCESS ) cout << "set  TRACE failed   status = "  <<  status << endl;

    params[0] = 100;
    status = opt->setOptimizerParameters( MAX_FUNCTION_EVALUATIONS, params );
    if( status != SUCCESS ) cout << "set  MAX_FUNCTION_EVALUATIONS failed   status = "  <<  status << endl;

    params[0] =  100;
    params[1] = -100;
    status = opt->setOptimizerParameters( INITIAL_VALUES,  params );
    if( status != SUCCESS ) cout << "set INITIAL_VALUES failed status = "  <<  status << endl;

    params[0] = .0001;
    status = opt->setOptimizerParameters( GRADIENT_CONVERGENCE_TOLERANCE, params );
    if( status != SUCCESS ) cout << "set  GRADIENT_CONVERGENCE_TOLERANCE failed   status = "  <<  status << endl;

    params[0] = 1.0;
    status = opt->setOptimizerParameters( DEFAULT_STEP_LENGTH, params );
    if( status != SUCCESS ) cout << "set  DEFAULT_STEP_LENGTH failed   status = "  <<  status << endl;

    params[0] = 0.9;
    status = opt->setOptimizerParameters( LINE_SEARCH_ACCURACY, params );
    if( status != SUCCESS ) cout << "set  LINE_SEARCH_ACCURACY failed   status = "  <<  status << endl;

    status = opt->setObjectiveFunction( costFunc );
    if( status != SUCCESS ) cout << "set  cost Function failed   status = "  <<  status << endl;

    status = opt->optimize( results );
    if( status != SUCCESS )  {
        cout << "Run Optimizer failed status = "   <<  status << endl;
    }else {
       for( i=0; i<PROBLEM_DIMENSION; i++ ) {
          printf(" results[%d] = %f \n",i,results[i]); 
       }
    }
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
