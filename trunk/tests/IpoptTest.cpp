/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKsimmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/BigMatrix.h"
#include "Optimizer.h"

#include <iostream>
using std::cout;
using std::endl;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;


static int  NUMBER_OF_PARAMETERS = 4; 
static int  NUMBER_OF_CONSTRAINTS = 2; 

/*
 * Adapted from Ipopt's hs071 example problem
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

class ProblemSystem : public OptimizerSystem {
public:


   int objectiveFunc(  const Vector &coefficients, const bool new_coefficients, Real& f ) const {
      const Real *x;

      x = &coefficients[0];

      f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, const bool new_coefficients, Vector &gradient ) const{
      const Real *x;

      x = &coefficients[0]; 

     gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
     gradient[1] = x[0] * x[3];
     gradient[2] = x[0] * x[3] + 1;
     gradient[3] = x[0] * (x[0] + x[1] + x[2]);

     return(0);

  }
  int constraintFunc( const Vector &coefficients, const bool new_coefficients, Vector &constraints)  const{
      const Real *x;

      x = &coefficients[0]; 
      constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
      constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

      return(0);
  }

  int constraintJacobian( const Vector& coefficients, const bool new_coefficients, Matrix& jac)  const{
//  int constraintJacobian( const Vector& coefficients, const bool new_coefficients, Vector& jac)  const{
      const Real *x;

      x = &coefficients[0]; 

      jac(0,0) = 2*x[0]; 
      jac(0,1) = 2*x[1];
      jac(0,2) = 2*x[2]; 
      jac(0,3) = 2*x[3]; 
      jac(1,0) = x[1]*x[2]*x[3]; 
      jac(1,1) = x[0]*x[2]*x[3];
      jac(1,2) = x[0]*x[1]*x[3]; 
      jac(1,3) = x[0]*x[1]*x[2]; 

      return(0);
  }

/*   ProblemSystem() : OptimizerSystem( NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS ) {} */

   ProblemSystem( const int numParams, const int numConstraints) :

         OptimizerSystem( numParams, numConstraints ) {
   }

};


int main() {

    Real f;
    int i;

    /* create the system to be optimized */
    ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS );

    Vector results(NUMBER_OF_PARAMETERS);
    Vector lower_bounds(NUMBER_OF_PARAMETERS);
    Vector upper_bounds(NUMBER_OF_PARAMETERS);


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


    int returnValue = 0; // assume success

  try {

    Optimizer opt( sys ); 


//    opt.setConvergenceTolerance( .0001 );
    opt.setConvergenceTolerance( 1e-3 );

    opt.setDiagnosticsLevel( 7 );

    /* compute  optimization */ 
    f = opt.optimize( results );
  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    returnValue = 1; // failure
    printf("IpoptTest.cpp: Caught exception \n");
  }

    printf("IpoptTest.cpp: f = %f params = ",f);
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");

    static const Real TOL = 1e-4;
    Real expected[] = { 1.00000000, 4.74299963, 3.82114998, 1.37940829 };
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       if( results[i] > expected[i]+TOL || results[i] < expected[i]-TOL) {
           printf(" IpoptTest.cpp: error results[%d] = %f  expected=%f \n",i,results[i], expected[i]);
           returnValue = 1;
       }
    }

    return( returnValue );
}
