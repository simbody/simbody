#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"

#include <iostream>
using std::cout;
using std::endl;


#define PROBLEM_DIMENSION 25
#define NUM_CORRECTIONS    5

class ProblemSystem : public SimTK::OptimizerSystem {

   int objectiveFunc(  int n, SimTK::Vector &coefficients, bool new_coefficients,  double *f  ) const  {
      double *x;
      int i;

      x = &coefficients[0];

//printf("objectiveFunction x = ",x[0],x[1],x[2]);
      *f = .25 *(x[0]-1.0)*(x[0]-1.0);
//   printf(" %f",x[0]);
      for(i=1;i<n;i++) {
         *f = *f + pow(x[i]-x[i-1]*x[i-1], 2.0);
//   printf(" %f",x[i]);
      }

//   printf(" \n");
      *f = 4.0* *f;
      return( 0 ); 
   }

   int gradientFunc(int n, SimTK::Vector &coefficients, bool new_coefficients,  SimTK::Vector &gradient ) const {
      double *x,t1,t2;
      int i;

      x = &coefficients[0]; 

      t1 = x[1]-(x[0]*x[0]);
      gradient[0] = 2.0*(x[0]-1.0)-16.0*x[0]*t1;
      for(i=1;i<n-1;i++) {
         t2=t1;
         t1=x[i+1]-(x[i]*x[i]);
         gradient[i]=8.0*t2-16.0*x[i]*t1;
      }
      gradient[n-1]=8.0*t1;
// printf("objectiveGradient x = %f %f %f  g = %f \n",x[0],x[1],x[2],gradient[0]);

    return(0);


   }

};

/* adapted from driver1.f of Lbfgsb.2.1.tar.gz  */
main() {

    double params[10],f;
    int bounds[PROBLEM_DIMENSION];
    int i;
    int n = PROBLEM_DIMENSION;

    SimTK::Vector results(PROBLEM_DIMENSION);
    SimTK::Vector lower_bounds(PROBLEM_DIMENSION);
    SimTK::Vector upper_bounds(PROBLEM_DIMENSION);
    ProblemSystem sys;

    cout << "LBFGSB driver1 test " << endl;
    sys.dimension = PROBLEM_DIMENSION;
    sys.numBounds = PROBLEM_DIMENSION;  // all coeffcients have bounds
    sys.lower_bounds = &lower_bounds[0];
    sys.upper_bounds = &upper_bounds[0];

    /* set initial conditions */
    for(i=0;i<n;i++) {
       results[i] = 3.0;
    }

    /* set bounds */
    for(i=0;i<n;i=i+2) {   // even numbered 
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 100.0;
       bounds[i] = i;
    }
    for(i=1;i<n;i=i+2) { // odd numbered
       lower_bounds[i] = -100.0;
       upper_bounds[i] = 100.0;
       bounds[i] = i;
    }


    try {
    SimTK::Optimizer opt( sys ); 

    params[0] = 100;
    opt.setOptimizerParameters( MAX_FUNCTION_EVALUATIONS, params );

    params[0] = .0001;
    opt.setOptimizerParameters( GRADIENT_CONVERGENCE_TOLERANCE, params );

    params[0] = 1.0;
    opt.setOptimizerParameters( DEFAULT_STEP_LENGTH, params );

    params[0] = 0.9;
    opt.setOptimizerParameters( LINE_SEARCH_ACCURACY, params );


    f = opt.optimize( results );

    }

    catch (SimTK::Exception::Base exp) {
        cout << "Caught exception :" << exp.getMessage() << endl;
    }

    printf("f = %f params = ",f);
    for( i=0; i<PROBLEM_DIMENSION; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");

}
