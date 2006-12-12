#ifndef _SimTK_LBFGSB_OPTIMIZER_H_
#define _SimTK_LBFGSB_OPTIMIZER_H_

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
#include "OptimizerRep.h"

namespace SimTK {


class LBFGSBOptimizer: public OptimizerRep {

     public:

     ~LBFGSBOptimizer() {

        delete [] nbd;

     }
     LBFGSBOptimizer(OptimizerSystem& sys); 
     unsigned int optParamStringToValue( char *parameter );
     void setOptimizerParameters(unsigned int parameter, double *values );
     void getOptimizerParameters(unsigned int parameter, double *values );
     double optimize(  Vector &results );

/* must implement get and set paramaeters and optimize() functions ?? optParamStringToValue ??*/
     private:
     int         numCorrections;
     int         MaxNumFuncEvals;
     double      GradientConvergenceTolerance;
     double      LineSearchAccuracy;
     double      DefaultStepLength;
     double      *gradient;
     int         iprint[2];
     int         *nbd;


};
} // namespace SimTK
#endif //_SimTK_LBFGSB_OPTIMIZER_H_

