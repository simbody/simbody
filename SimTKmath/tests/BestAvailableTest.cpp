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
 * Contributors: Chris Dembia                                                 *
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

/* This test doesn't actually run any optimizations. It simply creates an
 * Optimizer, using the same OptimizerSystem's (not the EXACT same) used in the
 * tests for each of the algorithms. NOTE: Since BestAvailable never selects
 * CFSQP, we do not test for CFSQP.
 */

/* See IpoptTest. 'BestAvailable' will select Ipopt because this sytem has
 * constraints.
 */
class IpoptSystem : public OptimizerSystem {
public:


    int objectiveFunc( const Vector &coefficients, bool new_coefficients,
            Real& f ) const {
        const Real *x;

        x = &coefficients[0];

        f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
        return( 0 ); 
    }

    int gradientFunc( const Vector &coefficients, bool new_coefficients,
            Vector &gradient ) const{
        const Real *x;

        x = &coefficients[0]; 

        gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        gradient[1] = x[0] * x[3];
        gradient[2] = x[0] * x[3] + 1;
        gradient[3] = x[0] * (x[0] + x[1] + x[2]);

        return(0);

    }

    int constraintFunc( const Vector &coefficients, bool new_coefficients,
            Vector &constraints ) const{
        const Real *x;

        x = &coefficients[0]; 
        constraints[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.0;
        constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

        return(0);
    }

    int constraintJacobian( const Vector& coefficients, bool new_coefficients,
            Matrix& jac ) const{
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


    IpoptSystem() : OptimizerSystem( 4 )
    {
        setNumEqualityConstraints( 1 );
        setNumInequalityConstraints( 1 );
    }

};


/* See LBFGSBTest.cpp. 'BestAvailable' will select LBFGSB because this system
 * does NOT have constraints and does have parameter limits.
 */
class LBFGSBSystem : public OptimizerSystem {

   public:

   LBFGSBSystem() : OptimizerSystem( 25 ) {}

   int objectiveFunc(   const Vector &coefficients, bool new_coefficients,  Real& f  ) const  {
      int i;

      const Real *x = &coefficients[0];

      f = .25 *(x[0]-1.0)*(x[0]-1.0);
      for(i=1;i<getNumParameters();i++) {
         f = f + pow(x[i]-x[i-1]*x[i-1], 2.0);
      }

      f = 4.0* f;
      return( 0 ); 
   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients,  Vector &gradient ) const {
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


/* See LBFGSTest.cpp. 'BestAvailable' will select LBFGS because this system
 * does NOT have constraints and does NOT have parameter limits.
 */
class LBFGSSystem : public OptimizerSystem {
   public:

   LBFGSSystem() : OptimizerSystem( 2 ) {}

   int objectiveFunc( const Vector &coefficients, bool new_coefficients,
           Real& f ) const {

      const Real x = coefficients[0];
      const Real y = coefficients[1];

      f = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y; 
    
      return(0);

   }

   int gradientFunc( const Vector &coefficients, bool new_coefficients,
           Vector &gradient ) const {

      const Real x = coefficients[0]; 
      const Real y = coefficients[1];  

      gradient[0] = 3*x + 2*y -2;
      gradient[1] = 2*x + 6*y +8; 

      return(0);

   }
};


int main() {

    int i;

    // IPOPT
    // -----
    /* create the system to be optimized */
    IpoptSystem sysIPOPT;
    Vector lowerBoundsIPOPT(sysIPOPT.getNumParameters());
    Vector upperBoundsIPOPT(sysIPOPT.getNumParameters());
    for(i=0;i<sysIPOPT.getNumParameters();i++) {   
        lowerBoundsIPOPT[i] = 1.0;
        upperBoundsIPOPT[i] = 5.0;
    }
    sysIPOPT.setParameterLimits(lowerBoundsIPOPT, upperBoundsIPOPT);


    // LBFGSB
    // ------
    LBFGSBSystem sysLBFGSB;
    Vector lowerBoundsLBFGSB(sysLBFGSB.getNumParameters());
    Vector upperBoundsLBFGSB(sysLBFGSB.getNumParameters());
    for(i=0;i<sysLBFGSB.getNumParameters();i=i+2) {   // even numbered 
       lowerBoundsLBFGSB[i] = 1.0;
       upperBoundsLBFGSB[i] = 100.0;
    }
    for(i=1;i<sysLBFGSB.getNumParameters();i=i+2) { // odd numbered
       lowerBoundsLBFGSB[i] = -100.0;
       upperBoundsLBFGSB[i] = 100.0;
    }
    sysLBFGSB.setParameterLimits( lowerBoundsLBFGSB, upperBoundsLBFGSB );


    // LBFSGB
    // ------
    LBFGSSystem sysLBFGS;


    int returnValue = 0; // assume success

    try {

        // Optimizer will choose the BestAvailable algorithm.
        Optimizer optIPOPT( sysIPOPT ); 
        if (optIPOPT.getAlgorithm() != SimTK::InteriorPoint)
        {   returnValue = 1; }

        Optimizer optLBFGSB( sysLBFGSB ); 
        if (optLBFGSB.getAlgorithm() != SimTK::LBFGSB)
        {   returnValue = 1; }

        Optimizer optLBFGS( sysLBFGS ); 
        if (optLBFGS.getAlgorithm() != SimTK::LBFGS)
        {   returnValue = 1; }

    }
    catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        returnValue = 1; // failure
        printf("BestAvailableTest.cpp: Caught exception \n");
    }

    return( returnValue );
}
