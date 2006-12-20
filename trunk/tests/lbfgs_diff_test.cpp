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

#include "Simmath.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"

#include <iostream>
using std::cout;
using std::endl;
using SimTK::Vector;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;

/* adapted from itkLBFGSOptimizerTest.cxx */

const static int  NUMBER_OF_PARAMETERS = 2;

class ProblemSystem : public OptimizerSystem {
   public:

   ProblemSystem( int numParameters) : OptimizerSystem( numParameters){}

   int objectiveFunc(  const Vector &coefficients, const bool new_coefficients, Real& f ) const {

      const Real x = coefficients[0];
      const Real y = coefficients[1];

      f = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y; 
    
      return(0);

   }

   int gradientFunc(  const Vector &coefficients, const bool new_coefficients, Vector &gradient )const {

      const Real x = coefficients[0]; 
      const Real y = coefficients[1];  

      gradient[0] = 3*x + 2*y -2;
      gradient[1] = 2*x + 6*y +8; 

      return(0);

   }
};

main() {

    Real params[10];
    int i;

    ProblemSystem sys(NUMBER_OF_PARAMETERS);

    Vector results(NUMBER_OF_PARAMETERS);

    Optimizer opt( sys ); 


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

    opt.useNumericalGradient( true );

    results[0] =  100;
    results[1] = -100;
    
    opt.optimize( results );


    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" results[%d] = %f \n",i,results[i]); 
    }

    static const Real TOL = 1e-4;
    Real expected[] = { 2.0, -2.0 };
    bool fail = false;
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       if( results[i] > expected[i]+TOL || results[i] < expected[i]-TOL) {
           printf(" lbfgs_diff_test error results[%d] = %f  expected=%f \n",i,results[i], expected[i]); 
           fail = true;
       }
    }

    if( fail ) 
       exit(1);
    else
       exit(0);

}
