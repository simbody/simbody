#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"

#include <iostream>
using std::cout;
using std::endl;


static int  NUMBER_OF_PARAMETERS = 4; 
static int  NUMBER_OF_CONSTRAINTS = 2; 

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

class ProblemSystem : public SimTK::OptimizerSystem {
public:


   int objectiveFunc(  Vector &coefficients, bool new_coefficients, Real* f ) const {
      double *x;
      int i;

      x = &coefficients[0];

      *f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( 0 ); 
   }

   int gradientFunc( Vector &coefficients, bool new_coefficients, Vector &gradient ) const{
      double *x;

      x = &coefficients[0]; 

     gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
     gradient[1] = x[0] * x[3];
     gradient[2] = x[0] * x[3] + 1;
     gradient[3] = x[0] * (x[0] + x[1] + x[2]);

     return(0);

  }
  int constraintFunc( Vector &coefficients, bool new_coefficients, Vector &constraints)  const{
      double *x;

      x = &coefficients[0]; 
      constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
      constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

      return(0);
  }

  int constraintJacobian( Vector& coefficients, bool new_coefficients, Vector& jac)  const{
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

      return(0);
  }

/*   ProblemSystem() : OptimizerSystem( NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS ) {} */

   ProblemSystem( int numParams, int numConstraints) :

         SimTK::OptimizerSystem( numParams, numConstraints ) {
   }

};

} // namespace SimTK

main() {

    double params[10],f;
    int i;

    /* create the system to be optimized */
    SimTK::ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS );

    SimTK::Vector results(NUMBER_OF_PARAMETERS);
    SimTK::Vector lower_bounds(NUMBER_OF_PARAMETERS);
    SimTK::Vector upper_bounds(NUMBER_OF_PARAMETERS);


    sys.setNumEqualityConstraints( 1 );

    /* set initial conditions */
    results[0] = 1.0;
    results[1] = 5.0;
    results[2] = 5.0;
    results[3] = 1.0;

    /* set bounds */
    for(i=0;i<NUMBER_OF_PARAMETERS;i++) {   
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 5.0;
    }

    sys.setParameterLimits( lower_bounds, upper_bounds );

    SimTK::Optimizer opt( sys ); 


    params[0] = 100;
    opt.setOptimizerParameters( MAX_FUNCTION_EVALUATIONS, params );

    params[0] = .0001;
    opt.setOptimizerParameters( GRADIENT_CONVERGENCE_TOLERANCE, params );

    params[0] = 1.0;
    opt.setOptimizerParameters( DEFAULT_STEP_LENGTH, params );

    params[0] = 0.9;
    opt.setOptimizerParameters( LINE_SEARCH_ACCURACY, params );

    /* compute  optimization */ 
    f = opt.optimize( results );

    printf("f = %f params = ",f);
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");


}
