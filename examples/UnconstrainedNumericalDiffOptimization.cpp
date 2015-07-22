/* -------------------------------------------------------------------------- *
 *    Simbody(tm) Example: Unconstrained Optimization w/Numerical Gradient    *
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
};

int main() {
    ProblemSystem sys(NUMBER_OF_PARAMETERS);

    Vector results(NUMBER_OF_PARAMETERS);

    Real f = NaN;
    try {
       Optimizer opt( sys );

       opt.setConvergenceTolerance( .0001 );
       opt.useNumericalGradient( true );

       results[0] =  100;
       results[1] = -100;

       f = opt.optimize( results );
    }
    catch(const std::exception& e) {
       cout << "ParameterConstrainedOptimization.cpp Caught exception :"  << endl;
       cout << e.what() << endl;
    }

    printf(" Optimal solution: f = %f   parameters = %f %f \n",f,results[0],results[1]);

    return 0;
}
