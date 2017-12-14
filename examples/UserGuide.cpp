/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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


   int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const override {
      const Real *x;

      x = &coefficients[0];

      f = (x[0] - 5.0)*(x[0] - 5.0) + (x[1] - 1.0)*(x[1] - 1.0);
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients, Vector &gradient ) const override{
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
  int constraintFunc( const Vector &coefficients, bool new_coefficients, Vector &constraints)  const override{
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
  int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const override{
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
