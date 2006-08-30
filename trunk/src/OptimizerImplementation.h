#ifndef _SimTK_OPTIMIZER_IMPLEMENTATION_H_
#define _SimTK_OPTIMIZER_IMPLEMENTATION_H_

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
#include "OptimizerInterface.h"
#include "ObjectiveFunction.h"

namespace SimTK {

class OptimizerImplementation : public OptimizerInterface {
    public:
     OptimizerImplementation( int );
     OptimizerImplementation();
     unsigned int optParamStringToValue( char *parameter );
     void setOptimizerParameters(unsigned  int parameter, double *values); 
     void getOptimizerParameters(unsigned int parameter, double *values); 
     void setObjectiveFunction( void (*func)(int,double*,double*,double*,void*)) ;
     void setObjectiveFunction( SimTK::ObjectiveFunction *objFunc);

     void optimize( double *results ); 

     void optimize(  SimTK::Vector &results ); 


      ~OptimizerImplementation(){
          if(dimension > 0 ) {
             delete [] work;
             delete [] diag;
             delete gradient;
         }
     }

     int         dimension;  // dimension of the problem
     int         numCorrections;
     int         Trace;
     int         Algorithm;
     int         MaxNumFuncEvals;
     double     *work;
     double     *diag;
     double      GradientConvergenceTolerance;
     double      LineSearchAccuracy;
     double      DefaultStepLength;
     void        (*costFunction)(int, double*, double*, double*, void*);
     int         iprint[2]; 
     double      xtol[1]; 
     int         diagco[1];
     void       *user_data;
     SimTK::Vector *gradient;
     SimTK::ObjectiveFunction *objFunc;


};
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_IMPLEMENTATION_H_
