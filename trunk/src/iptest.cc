#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"
#include "OptimizationProblem.h"

#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {
/*
 * Problem hs071 looks like this
 *
 *     min   x1*x4*(x1 + x2 + x3)  +  x3
 *     s.t.  x1*x2*x3*x4                   >=  25
 *           x1**2 + x2**2 + x3**2 + x4**2  =  40
 *           1 <=  x1,x2,x3,x4  <= 5
 *
 *     Starting point:
 *        x = (1, 5, 5, 1)
 *
 *     Optimal solution:
 *        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
 *
 */
#define PROBLEM_DIMENSION 4 

class ProblemStatement : public SimTK::OptimizationProblem {

   double objectiveFunction(  int n, bool new_coefficients, Vector &coefficients, void* user_data ) {
      double *x,f;
      int i;

      x = &coefficients[0];

//printf("objectiveFunction x = ",x[0],x[1],x[2]);
      f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( f ); 
   }

   void objectiveGradient(int n,   bool new_coefficients, Vector &coefficients, Vector &gradient, void *user_data ){
      double *x;

      x = &coefficients[0]; 

     gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
     gradient[1] = x[0] * x[3];
     gradient[2] = x[0] * x[3] + 1;
     gradient[3] = x[0] * (x[0] + x[1] + x[2]);

  }
  void computeConstraints(int n, int m, bool new_coefficients, Vector & coefficients, Vector & constraints, void * user_data) {
      double *x;

      x = &coefficients[0]; 
          constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
          constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

           return;
  }

  void computeConstraintsJacobian(int n, int m, bool new_coefficients, Vector& coefficients, Vector& jac, void * user_data) {
      double *x;

      x = &coefficients[0]; 

    jac[0] = 2*x[0]; // 1,0
    jac[1] = 2*x[1]; // 1,1
    jac[2] = 2*x[2]; // 1,2
    jac[3] = 2*x[3]; // 1,3
    jac[4] = x[1]*x[2]*x[3]; // 0,0
    jac[5] = x[0]*x[2]*x[3]; // 0,1
    jac[6] = x[0]*x[1]*x[3]; // 0,2
    jac[7] = x[0]*x[1]*x[2]; // 0,3


          return;
  }

};
}

main() {

    double params[10],f;
    int i;
    int n = PROBLEM_DIMENSION;

    SimTK::Vector results(PROBLEM_DIMENSION);
    SimTK::Vector lower_bounds(PROBLEM_DIMENSION);
    SimTK::Vector upper_bounds(PROBLEM_DIMENSION);
    SimTK::ProblemStatement study;

    study.dimension = PROBLEM_DIMENSION;
    study.numBounds = PROBLEM_DIMENSION;  // all coeffcients have bounds
    study.lower_bounds = &lower_bounds[0];
    study.upper_bounds = &upper_bounds[0];
    study.numConstraints = 2;
    study.numEqualityConstraints = 1;

    /* set initial conditions */
    results[0] = 1.0;
    results[1] = 5.0;
    results[2] = 5.0;
    results[3] = 1.0;

    /* set bounds */
    for(i=0;i<n;i++) {   
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 5.0;
    }


    try {
    SimTK::Optimizer opt( study ); 

    params[0] = 0;
    opt.setOptimizerParameters( TRACE, params );

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
