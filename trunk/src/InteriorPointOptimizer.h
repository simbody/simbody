#ifndef _SimTK_INTERIOR_POINT_OPTIMIZER_H_
#define _SimTK_INTERIOR_POINT_OPTIMIZER_H_

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
//#include <iostream>
#include "Simmath.h"
#include "OptimizerInterface.h"
#include "OptimizationProblem.h"
#include "ipopt/IpStdCInterface.h"
#include "ipopt/IpReturnCodes.h"

namespace SimTK {


class InteriorPointOptimizer: public OptimizerInterface {

     public:

     ~InteriorPointOptimizer() {
        printf("InteriorPointOptimizer destructor \n");


     }
     InteriorPointOptimizer(OptimizationProblem& p); 
     InteriorPointOptimizer(); 
     unsigned int InteriorPointOptimizer::optParamStringToValue( char *parameter );
     void setOptimizerParameters(unsigned int parameter, double *values );
     void getOptimizerParameters(unsigned int parameter, double *values );
     double optimize( double *results );
     double optimize(  SimTK::Vector &results );
/*
     Bool f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

     Bool grad_f(Index n, Number* x, Bool new_x,
            Number* grad_func, UserDataPtr user_data);

     Bool g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

     Bool jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

     Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);
*/
     
/* must implement get and set paramaeters and optimize() functions ?? optParamStringToValue ??*/
     private:
     OptimizationProblem *pop;
     int         n;
     int         m;
     int         numCorrections;
     int         MaxNumFuncEvals;
     double      GradientConvergenceTolerance;
     double      LineSearchAccuracy;
     double      DefaultStepLength;
     Number      *mult_x_L;
     Number      *mult_x_U;
     Number      *g_L;
     Number      *g_U;
     IpoptProblem nlp;
     SimTK::Vector *gradient;


};
} // namespace SimTK
#endif //_SimTK_INTERIOR_POINT_OPTIMIZER_H_

