/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKmath.h"

#include <iostream>
using std::cout;
using std::endl;
using namespace SimTK;

/* adapted from itkLBFGSOptimizerTest.cxx */

const static int  NUMBER_OF_PARAMETERS = 2;

class ProblemSystem : public OptimizerSystem {
   public:

   ProblemSystem( int numParameters) : OptimizerSystem( numParameters){}

   int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const {

      const Real x = coefficients[0];
      const Real y = coefficients[1];

      f = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y; 
    
      return(0);

   }

   int gradientFunc(  const Vector &coefficients, bool new_coefficients, Vector &gradient )const {

      const Real x = coefficients[0]; 
      const Real y = coefficients[1];  

      gradient[0] = 3*x + 2*y -2;
      gradient[1] = 2*x + 6*y +8; 

      return(0);

   }
};

static const Real expected[] = { 2.0, -2.0 };


static bool equalToTol(Real v1, Real v2, Real tol) {
    const Real scale = std::max(std::max(std::abs(v1), std::abs(v2)), Real(1));
    return std::abs(v1-v2) < scale*tol;
}

int main() {

    int i;

    ProblemSystem sys(NUMBER_OF_PARAMETERS);

    Vector results(NUMBER_OF_PARAMETERS);

      int returnValue = 0; // assume success
  try {

    
    Optimizer opt( sys ); 

    opt.setConvergenceTolerance( .0001 );
    results[0] =  100;
    results[1] = -100;
    
    opt.setAdvancedRealOption( "xtol", 1e-6 );

    opt.optimize( results );

  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    returnValue = 1; // failure
  }



    static const Real TOL = 1e-4;

    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       if(!equalToTol(results[i], expected[i], TOL)) {
           printf(" LBFGSTest.cpp:  error results[%d] = %f  expected=%f \n",
                  i,results[i], expected[i]); 
           returnValue = 1;
       }
    }

    if( returnValue == 0 ) {
        printf("LBFGSTest.cpp results = ");
        for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
             printf( "%f ", results[i]);
        }
        printf("\n");
    }

    return( returnValue );

  
}
