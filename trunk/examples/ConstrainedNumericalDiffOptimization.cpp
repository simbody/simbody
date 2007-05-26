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

#include <iostream>
#include "Simmath.h"
#include "Optimizer.h"

using namespace SimTK;

static int  NUMBER_OF_PARAMETERS = 4; 
static int  NUMBER_OF_EQUALITY_CONSTRAINTS = 1; 
static int  NUMBER_OF_INEQUALITY_CONSTRAINTS = 1; 

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

class ProblemSystem : public OptimizerSystem {
public:


   int objectiveFunc(  const Vector &coefficients, const bool new_coefficients, Real& f ) const {
      const Real *x;
      int i;

      x = &coefficients[0];

      f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( 0 ); 
   }

  int constraintFunc( const Vector &coefficients, const bool new_coefficients, Vector &constraints)  const{
      const Real *x;

      x = &coefficients[0]; 
      constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
      constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

      return(0);
  }

    ProblemSystem( const int numParams, const int numEqualityConstraints, const int numInequalityConstraints ) :
        OptimizerSystem( numParams ) 
    {
        setNumEqualityConstraints( numEqualityConstraints );
        setNumInequalityConstraints( numInequalityConstraints );
    }

};


main() {

    Real f;
    int i;

    /* create the system to be optimized */
    ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_EQUALITY_CONSTRAINTS, NUMBER_OF_INEQUALITY_CONSTRAINTS);

    Vector results(NUMBER_OF_PARAMETERS);
    Vector lower_bounds(NUMBER_OF_PARAMETERS);
    Vector upper_bounds(NUMBER_OF_PARAMETERS);

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

    try {
       Optimizer opt( sys ); 

       opt.setConvergenceTolerance( .0001 );
       opt.useNumericalGradient( true );
       opt.useNumericalJacobian( true );

       /* compute  optimization */ 
       f = opt.optimize( results );
    }
    catch (const std::exception& e) {
       std::cout << "ConstrainedNumericalDiffOptimization.cpp: Caught exception" << std::endl;
       std::cout << e.what() << std::endl;
    }

    printf("Optimal Solution: f = %f   parameters = %f %f %f %f \n",f, results[0],results[1],results[2],results[3]);

}
