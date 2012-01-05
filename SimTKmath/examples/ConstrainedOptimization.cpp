
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

#include "SimTKmath.h"
#include <iostream>

using namespace SimTK;

static int  NumberOfParameters = 4; 
static int  NumberOfEqualityConstraints = 1; 
static int  NumberOfInequalityConstraints = 1; 

/*
 * This example was adapted from IPOPT's hs071 example 
 *
 *   Problem statement:
 *
 *     minimize:   x1*x4*(x1 + x2 + x3)  +  x3
 *
 *     s.t.  x1*x2*x3*x4                   >=  25    inequality constraint 
 *           x1**2 + x2**2 + x3**2 + x4**2 - 40.0 = 0.0  equality constraint
 *
 *           1 <=  x1,x2,x3,x4  <= 5    each parameter has a lower limit of 1.0 and an upper limit of 5.0
 *
 *     Starting point:
 *        x = (1, 5, 5, 1)   will be used for the initial conditions
 *
 *     Optimal solution:
 *        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
 *
 */

class ProblemSystem : public OptimizerSystem {
public:


   int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const {
      const Real *x;

      x = &coefficients[0];

      f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients, Vector &gradient ) const{
      const Real *x;

      x = &coefficients[0]; 

     gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
     gradient[1] = x[0] * x[3];
     gradient[2] = x[0] * x[3] + 1;
     gradient[3] = x[0] * (x[0] + x[1] + x[2]);

     return(0);

  }

  /* 
  ** Method to compute the value of the constraints.
  ** Equality constraints are first followed by the any inequality constraints
  */ 
  int constraintFunc( const Vector &coefficients, bool new_coefficients, Vector &constraints)  const{
      const Real *x;

      x = &coefficients[0]; 
      constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
      constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

      return(0);
  }


  /*
  ** Method to compute the jacobian of the constraints.
  **
  */
  int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const{
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

    ProblemSystem( const int numParams, const int numEqualityConstraints, const int numInequalityConstraints ) :
        OptimizerSystem( numParams ) 
    {
        setNumEqualityConstraints( numEqualityConstraints );
        setNumInequalityConstraints( numInequalityConstraints );
    }

};


int main() {
    /* create the system to be optimized */
    ProblemSystem sys(NumberOfParameters, NumberOfEqualityConstraints, NumberOfInequalityConstraints);

    Vector results(NumberOfParameters);
    Vector lower_bounds(NumberOfParameters);
    Vector upper_bounds(NumberOfParameters);

    /* set initial conditions */
    results[0] = 1.0;
    results[1] = 5.0;
    results[2] = 5.0;
    results[3] = 1.0;

    /* set bounds */
    for(int i=0;i<NumberOfParameters;i++) {   
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 5.0;
    }

    sys.setParameterLimits( lower_bounds, upper_bounds );

    Real f = NaN;
    try {
        Optimizer opt( sys ); 

        opt.setConvergenceTolerance( .0001 );

        /* compute  optimization */ 
        f = opt.optimize( results );
    }
    catch (const std::exception& e) {
        std::cout << "ConstrainedOptimization.cpp Caught exception:" << std::endl;
        std::cout << e.what() << std::endl;
    }


    printf("Optimal Solution: f = %f   parameters = %f %f %f %f \n",f,results[0],results[1],results[2],results[3]);

    return 0;
}
