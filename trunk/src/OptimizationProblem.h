#ifndef _SimTK_OPTIMIZATION_PROBLEM_H
#define _SimTK_OPTIMIZATION_PROBLEM_H

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

namespace SimTK {

class OptimizationProblem { 

public:
    OptimizationProblem() { 
       printf("default OptimizationProblem constructor \n");
       dimension = 1; 
       numConstraints = 0;
       numEqualityConstraints = 0;
    }
    OptimizationProblem(int n) { 
       printf("OptimizationProblem constructor n=%d \n",n);
       dimension = n; 
       numConstraints = 0;
       numEqualityConstraints = 0;
    }

    virtual ~OptimizationProblem() {};

  /* this method must be overloaded by derived class */

  virtual double objectiveFunction( int n, 
                                    SimTK::Vector& coefficients, 
                                    void* user_data  ) = 0;

  virtual void objectiveGradient(int n,  
                                 SimTK::Vector &coefficients,   
                                 SimTK::Vector &gradient, 
                                 void* user_data) {};

  virtual void computeConstraints(int n,  
                                  int m,
                                  SimTK::Vector & params,   
                                  SimTK::Vector & constraints, 
                                  void * user_data) {
           return;
  }

  virtual void computetConstraintJacobian(int n,  
                                          int m,
                                          SimTK::Vector & params,   
                                          SimTK::Vector & jac, 
                                          void * user_data) {
          return;
  }

  virtual double objectiveAndGradient( int n,
                                       SimTK::Vector &coefficients, 
                                       SimTK::Vector &gradient,
                                       void *user_data) {

     this->objectiveGradient( n, coefficients, gradient, user_data );
     return( this->objectiveFunction(n, coefficients, user_data ) );

  }
  virtual bool computeHessian( int n,
                               int m, 
                               SimTK::Vector &coefficients, 
                               SimTK::Vector &gradient,
                               void *user_data) { 
          return 1;
 }



   
      int dimension;    // number of indepent varaibles
      int numConstraints;    // number of constraints
      int numEqualityConstraints;
      void *user_data;


}; // end class OptimizationProblem
}

#endif //_SimTK_OPTIMIZATION_PROBLEM_H



