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
using SimTK::Vector;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;


const static int NUMBER_OF_PARAMETERS = 25;

class ProblemSystem : public OptimizerSystem {

   public:

   ProblemSystem( const int numParameters ) : OptimizerSystem( numParameters ) {}

   int objectiveFunc(   const Vector &coefficients, bool new_coefficients,  Real& f  ) const override  {
      int i;

      const Real *x = &coefficients[0];

      f = .25 *(x[0]-1.0)*(x[0]-1.0);
      for(i=1;i<getNumParameters();i++) {
         f = f + pow(x[i]-x[i-1]*x[i-1], 2.0);
      }

      f = 4.0* f;
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients,  Vector &gradient ) const override {
      const Real *x;
      Real t1,t2;
      int i;

      x = &coefficients[0]; 

      t1 = x[1]-(x[0]*x[0]);
      gradient[0] = 2.0*(x[0]-1.0)-16.0*x[0]*t1;
      for(i=1;i<getNumParameters()-1;i++) {
         t2=t1;
         t1=x[i+1]-(x[i]*x[i]);
         gradient[i]=8.0*t2-16.0*x[i]*t1;
      }
      gradient[getNumParameters()-1]=8.0*t1;

    return(0);

   }

};

static bool equalToTol(Real v1, Real v2, Real tol) {
    const Real scale = std::max(std::max(std::abs(v1), std::abs(v2)), Real(1));
    return std::abs(v1-v2) < scale*tol;
} 

/* adapted from driver1.f of Lbfgsb.2.1.tar.gz  */
int main() {

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

  int returnValue = 0; // assume success

  try {
    Optimizer opt( sys ); 

    opt.setConvergenceTolerance( .0001 );
    opt.setDiagnosticsLevel( 5 );
    opt.setAdvancedRealOption( "factr", 1000);
    f = opt.optimize( results );

  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    returnValue = 1; // failure
  }



    printf("LBFGSBTest.cpp: f = %f params = ",f);
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       printf(" %f",results[i]); 
    }
    printf("\n");

    static const Real TOL = 1e-3;
    Real expected[] = { 1.000000, 0.999998, 1.000000, 1.000001, 1.000003, 
                        1.000006, 1.000007, 1.000012, 1.000022, 1.000040, 
                        1.000081, 1.000161, 1.000325, 1.000650, 1.001302, 
                        1.002603, 1.005214, 1.010450, 1.021013, 1.042466, 
                        1.086736, 1.180997, 1.394759, 1.945352, 3.784388 };
    for( i=0; i<NUMBER_OF_PARAMETERS; i++ ) {
       if(!equalToTol(results[i], expected[i], TOL)) {
           printf(" LBFGSBTest.cpp: error results[%d] = %f  expected=%f \n",i,results[i], expected[i]);
           returnValue = 1;
       }
    }

    return( returnValue );


}
