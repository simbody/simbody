
#ifndef _SimTK_OPTIMIZER_H
#define _SimTK_OPTIMIZER_H

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



#include "OptimizerInterface.h"
#include "OptimizationProblem.h"
#include "ObjectiveFunction.h"

namespace SimTK {
/*
** Class for API interface to Simmath's optimizers.
** The OptimizationProblem class describes the optimization by
** specifying the objective function and constraints. OptimizerFactory()
** instantiates the coorect optimizer based on the objective function 
** and constraints specified in the OptimizationProblem object. 
** If the user calls the Optimizer constructor and 
** supplies the algorithm argument the OptimizerFactory() will ignore the 
** will create instatiate the Optimizer asked for.
**  
*/

class Optimizer :  public OptimizerInterface {

   public:
    Optimizer( OptimizationProblem& op);
    Optimizer( OptimizationProblem& op, OptimizerAlgorithm opt_algo);
    Optimizer( int n, int nConstraints, int nEqualConstraints, int nBounds);

    ~Optimizer() {
       delete( (OptimizerInterface *)data );
    }
    void setOptimizerParameters(unsigned int param, double *values); 
    void getOptimizerParameters(unsigned int param, double *values);

    double optimize(SimTK::Vector&) ;

    private:
     void *data;
     OptimizerInterface *OptimizerFactory(OptimizationProblem& );
     OptimizerInterface *OptimizerFactory(OptimizationProblem&, OptimizerAlgorithm );
     OptimizerInterface *OptimizerFactory(int n, int nConstraints, int nEqualConstraints, int nBounds);

}; // Class Optimizer
} // namespace SimTK

#endif //_SimTK_OPTIMIZER_H

