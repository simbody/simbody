
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
#include "Optimizer.h"
#include "OptimizerImplementation.h"
#include "OptimizationProblem.h"
#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"

namespace SimTK {
   Optimizer::Optimizer(OptimizationProblem& problem) {
        OptimizerImplementation *optPtr =  OptimizerFactory(problem);
        data = (void *)optPtr;
   }
   Optimizer::Optimizer(OptimizationProblem& problem, OptimizerAlgorithm algorithm) {
        OptimizerImplementation *optPtr =  OptimizerFactory(problem, algorithm);
        data = (void *)optPtr;
   }

   void  Optimizer::setOptimizerParameters(unsigned int param, double *values) {

 //     ((OptimizerImplementation *)data)->setOptimizerParameters(param, values);
      return;
   }
   OptimizerImplementation *Optimizer::OptimizerFactory( OptimizationProblem& problem, OptimizerAlgorithm algorithm) {
 /*  
     if( algrothim == LBFGS) {
        return (OptimizerImplementation *) new LBFGSOptimizer( problem  );
     } else if( algrothim == LBFGSB) {
        return (OptimizerImplementation *) new LBFGSBOptimizer( problem  );
     } else {
        return (OptimizerImplementation *) new InteriorPointOptimizer( problem  );
     }
*/
  }
  OptimizerImplementation *Optimizer::OptimizerFactory( OptimizationProblem& problem) {
/*   
     if( problem.nconstraints > 0   {
        return (OptimizerImplementation *) new InteriorPointOptimizer( problem  );
*/
     if( problem.numBounds > 0 ) {
        return (OptimizerImplementation *) new LBFGSBOptimizer( problem  );
     } else {
        return (OptimizerImplementation *) new LBFGSOptimizer( problem  );
     }

  }

   void Optimizer::getOptimizerParameters(unsigned int param, double *values) {

  //    ((OptimizerImplementation *)data)->getOptimizerParameters(param, values);
      return;
   }

   void Optimizer::optimize(SimTK::Vector   &results) {
      ((OptimizerImplementation *)data)->optimize(results);
       return;
   }

} // namespace SimTK
