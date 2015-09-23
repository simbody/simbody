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
using SimTK::Matrix;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;


static int  NUMBER_OF_PARAMETERS = 4;
static int  NUMBER_OF_EQUALITY_CONSTRAINTS = 1;
static int  NUMBER_OF_INEQUALITY_CONSTRAINTS = 1;

/*
 * Adapted from Ipopt's hs071 example
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


   int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const override {
      const Real *x;

      x = &coefficients[0];

      f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
      return( 0 );
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients, Vector &gradient ) const override{
      const Real *x;

      x = &coefficients[0];

     gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
     gradient[1] = x[0] * x[3];
     gradient[2] = x[0] * x[3] + 1;
     gradient[3] = x[0] * (x[0] + x[1] + x[2]);

     return(0);

  }
  int constraintFunc( const Vector &coefficients, bool new_coefficients, Vector &constraints)  const override{
      const Real *x;

      x = &coefficients[0];
      constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
      constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

      return(0);
  }

  int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const override{
      const Real *x;

      x = &coefficients[0];

      jac[0][0] = 2*x[0];
      jac[0][1] = 2*x[1];
      jac[0][2] = 2*x[2];
      jac[0][3] = 2*x[3];
      jac[1][0] = x[1]*x[2]*x[3];
      jac[1][1] = x[0]*x[2]*x[3];
      jac[1][2] = x[0]*x[1]*x[3];
      jac[1][3] = x[0]*x[1]*x[2];


      return(0);
  }

/*   ProblemSystem() : OptimizerSystem( NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS ) {} */

    ProblemSystem( const int numParams, const int numEqualityConstraints, const int numInequalityConstraints ) :
        OptimizerSystem( numParams )
    {
        setNumEqualityConstraints( numEqualityConstraints );
        setNumInequalityConstraints( numInequalityConstraints );
   }

};


int main() {

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

    int returnValue = 0; // assume success
  try {

    Optimizer opt( sys );

    opt.setConvergenceTolerance( 1e-4 );

    opt.useNumericalGradient( true );
    opt.useNumericalJacobian( true );
    opt.setDiagnosticsLevel( 5 );

    /* compute  optimization */
    f = opt.optimize( results );
  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    returnValue = 1; // failure
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
