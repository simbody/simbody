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

const static int NUMBER_OF_PARAMETERS = 25;

class ProblemSystem : public OptimizerSystem {

   public:

   ProblemSystem( const int numParameters ) : OptimizerSystem( numParameters ) {}

   int objectiveFunc(   const Vector &coefficients, const bool new_coefficients,  Real& f  ) const  {
      int i;

      const Real *x = &coefficients[0];

      f = .25 *(x[0]-1.0)*(x[0]-1.0);
      for(i=1;i<numParameters;i++) {
         f = f + pow(x[i]-x[i-1]*x[i-1], 2.0);
      }

      f = 4.0* f;
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, const bool new_coefficients,  Vector &gradient ) const {
      const Real *x;
      Real t1,t2;
      int i;

      x = &coefficients[0]; 

      t1 = x[1]-(x[0]*x[0]);
      gradient[0] = 2.0*(x[0]-1.0)-16.0*x[0]*t1;
      for(i=1;i<numParameters-1;i++) {
         t2=t1;
         t1=x[i+1]-(x[i]*x[i]);
         gradient[i]=8.0*t2-16.0*x[i]*t1;
      }
      gradient[numParameters-1]=8.0*t1;

    return(0);

   }

};

/* adapted from driver1.f of Lbfgsb.2.1.tar.gz  */
main() {

    Real f;
    int i;
    int n = NUMBER_OF_PARAMETERS;

    Vector results(NUMBER_OF_PARAMETERS);
    Vector lower_bounds(NUMBER_OF_PARAMETERS);
    Vector upper_bounds(NUMBER_OF_PARAMETERS);

    ProblemSystem sys(NUMBER_OF_PARAMETERS);


    /* set initial conditions */
    for(i=0;i<n;i++) {
       results[i] = 3.0;
    }

    /* set bounds */
    for(i=0;i<n;i=i+2) {   // even numbered 
       lower_bounds[i] = 1.0;
       upper_bounds[i] = 100.0;
    }
    for(i=1;i<n;i=i+2) { // odd numbered
       lower_bounds[i] = -100.0;
       upper_bounds[i] = 100.0;
    }

    sys.setParameterLimits( lower_bounds, upper_bounds );

    try {
       Optimizer opt( sys ); 

       opt.setConvergenceTolerance( .0001 );

       f = opt.optimize( results );

    }

    catch (SimTK::Exception::Base e) {
        std::cout << "ParameterConstrainedOptimization.cpp Caught exception :" <<  std::endl;
        std::cout << e.what() << std::endl;
    }

    printf("f = %f params = ",f);
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");

}
