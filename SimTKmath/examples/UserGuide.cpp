
/* Portions copyright (c) 2007 Stanford University and Jack Middleton.
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

static int  NUMBER_OF_PARAMETERS = 2; 
static int  NUMBER_OF_EQUALITY_CONSTRAINTS = 0; 
static int  NUMBER_OF_INEQUALITY_CONSTRAINTS = 2; 

/*
 *
 *   Problem statement:
 *
 *     minimize:   (x1+5)^2 + (x2 - 4)^2 
 *
 *     s.t.   x1 - x2^2 <= 0    // inequality constraint
              x1 + x2 >= -2     //    inequality constraint 
 *
 *     Starting point:
 *        x = ( 5, 5)   will be used for the initial conditions
 *
 *     Optimal solution:
 *        x = (1.00000000 )
 *
 */

class ProblemSystem : public OptimizerSystem {
public:


   int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const {
      const Real *x;

      x = &coefficients[0];

      f = (x[0] - 5.0)*(x[0] - 5.0) + (x[1] - 1.0)*(x[1] - 1.0);
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients, Vector &gradient ) const{
      const Real *x;

      x = &coefficients[0]; 

     gradient[0] = 2.0*(x[0] - 5.0);
     gradient[1] = 2.0*(x[1] - 1.0);

     return(0);

  }

  /* 
  ** Method to compute the value of the constraints.
  ** Equality constraints are first followed by the any inequality constraints
  */ 
  int constraintFunc( const Vector &coefficients, bool new_coefficients, Vector &constraints)  const{
      const Real *x;

      x = &coefficients[0]; 
      constraints[0] = x[0] - x[1]*x[1];
      constraints[1] = x[1] - x[0] + 2.0;

      return(0);
  }


  /*
  ** Method to compute the jacobian of the constraints.
  **
  */
  int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const{
      const Real *x;

      x = &coefficients[0]; 
      jac(0,0) =  1.0;
      jac(0,1) = -2.0*x[1];
      jac(1,0) = -1.0;
      jac(1,1) =  1.0;


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
    ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_EQUALITY_CONSTRAINTS, NUMBER_OF_INEQUALITY_CONSTRAINTS);

    Vector results(NUMBER_OF_PARAMETERS);

    /* set initial conditions */
    results[0] = 5.0;
    results[1] = 5.0;

    Real f = NaN;

    try {
        Optimizer opt( sys ); 

        opt.setConvergenceTolerance( .0000001 );

        /* compute  optimization */ 
        f = opt.optimize( results );
    }
    catch (const std::exception& e) {
        std::cout << "ConstrainedOptimization.cpp Caught exception:" << std::endl;
        std::cout << e.what() << std::endl;
    }


    printf("Optimal Solution: f = %f   parameters = %f %f \n",f,results[0],results[1]);

    return 0;
}
