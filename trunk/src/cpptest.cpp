#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "optimizer.h"
#include "objectiveFunction.h"

#include <iostream>
using std::cout;
using std::endl;


#define PROBLEM_DIMENSION 2

class objfunc : public SimTK::objectiveFunction {

   double getValue(  SimTK::Vector_<SimTK::Real> &coefficients ) {
      double x, y;

      x = coefficients[0]; 
      y = coefficients[1];  

      return( 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y); 

   }

   void getGradient(  SimTK::Vector_<SimTK::Real> &coefficients, SimTK::Vector_<SimTK::Real> &gradient ){
      double x, y;

      x = coefficients[0]; 
      y = coefficients[1];  

      gradient[0] = 3*x + 2*y -2;
      gradient[1] = 2*x + 6*y +8; 
   }

};

/* adapted from itkLBFGSOptimizerTest.cxx */
main() {

    double params[10];
    int status,i;

    objfunc of;

    
    cout << "cpptest " << endl;
    SimTK::smOptimizer *opt = new SimTK::smOptimizer( PROBLEM_DIMENSION );

    params[0] = 0;
    status = opt->setOptimizerParameters( TRACE, params );
    if( status != SUCCESS ) cout << "set  TRACE failed   status = "  <<  status << endl;

    params[0] = 100;
    status = opt->setOptimizerParameters( MAX_FUNCTION_EVALUATIONS, params );
    if( status != SUCCESS ) cout << "set  MAX_FUNCTION_EVALUATIONS failed   status = "  <<  status << endl;


    params[0] = .0001;
    status = opt->setOptimizerParameters( GRADIENT_CONVERGENCE_TOLERANCE, params );
    if( status != SUCCESS ) cout << "set  GRADIENT_CONVERGENCE_TOLERANCE failed   status = "  <<  status << endl;

    params[0] = 1.0;
    status = opt->setOptimizerParameters( DEFAULT_STEP_LENGTH, params );
    if( status != SUCCESS ) cout << "set  DEFAULT_STEP_LENGTH failed   status = "  <<  status << endl;

    params[0] = 0.9;
    status = opt->setOptimizerParameters( LINE_SEARCH_ACCURACY, params );
    if( status != SUCCESS ) cout << "set  LINE_SEARCH_ACCURACY failed   status = "  <<  status << endl;

//    status = opt->setObjectiveFunction( costFunc );

    status = opt->setObjectiveFunction( &of );
    if( status != SUCCESS ) cout << "set  cost Function failed   status = "  <<  status << endl;

    SimTK::Vector_<SimTK::Real> results(2);
    results[0] =  100;
    results[1] = -100;
    
    status = opt->optimize( results );
    if( status != SUCCESS )  {
        cout << "Run Optimizer failed status = "   <<  status << endl;
    }else {
       for( i=0; i<PROBLEM_DIMENSION; i++ ) {
          printf(" results[%d] = %f \n",i,results[i]); 
       }
    }
}
