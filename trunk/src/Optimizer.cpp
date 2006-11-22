
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
#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"
#include "InteriorPointOptimizer.h"
#include "OptimizationProblem_C.h"

namespace SimTK {
   Optimizer::Optimizer(OptimizationProblem& problem) {
        OptimizerInterface *optPtr =  OptimizerFactory(problem);
        data = (void *)optPtr;
   }
   Optimizer::Optimizer(int dimension, int nConstraints, int nEqualConstraints, int nBounds) {
        OptimizationProblem_C problem( dimension, nConstraints, nEqualConstraints, nBounds );
        
        OptimizerInterface *optPtr =  OptimizerFactory(problem);
        data = (void *)optPtr;
   }
   Optimizer::Optimizer(OptimizationProblem& problem, OptimizerAlgorithm algorithm) {
        OptimizerInterface *optPtr =  OptimizerFactory(problem, algorithm);
        data = (void *)optPtr;
   }

   void  Optimizer::setOptimizerParameters(unsigned int param, double *values) {

      ((OptimizerInterface *)data)->setOptimizerParameters(param, values);
      return;
   }
   OptimizerInterface *Optimizer::OptimizerFactory( OptimizationProblem& problem, OptimizerAlgorithm algorithm) {
   
     if( algorithm == LBFGS) {
        return (OptimizerInterface *) new LBFGSOptimizer( problem  );
     } else if( algorithm == LBFGSB) {
        return (OptimizerInterface *) new LBFGSBOptimizer( problem  );
     } else {
        return (OptimizerInterface *) new InteriorPointOptimizer( problem  );
     }

  }
  OptimizerInterface *Optimizer::OptimizerFactory( OptimizationProblem& problem) {
   
     if( problem.numConstraints > 0)   {
        return (OptimizerInterface *) new InteriorPointOptimizer( problem  );
     }else if( problem.numBounds > 0 ) {
        return (OptimizerInterface *) new LBFGSBOptimizer( problem  );
     } else {
        return (OptimizerInterface *) new LBFGSOptimizer( problem  );
     }

  }

   void Optimizer::getOptimizerParameters(unsigned int param, double *values) {

      ((OptimizerInterface *)data)->getOptimizerParameters(param, values);
      return;
   }

   double Optimizer::optimize(SimTK::Vector   &results) {
      ((OptimizerInterface *)data)->optimize(results);
       return(0.0);
   }

} // namespace SimTK
