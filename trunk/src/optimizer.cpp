#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

namespace SimTK {
   smOptimizer::smOptimizer(int dimension) {
        optimizerImplementation *optPtr =  new optimizerImplementation(dimension);
        optPtr->dimension = dimension;
        data = (void *)optPtr;
   }

   smStatus  smOptimizer::setOptimizerParameters(unsigned int param, double *values) {

         return(((optimizerImplementation *)data)->setOptimizerParameters(param, values));
   }

   smStatus smOptimizer::getOptimizerParameters(unsigned int param, double *values) {

         return(((optimizerImplementation *)data)->getOptimizerParameters(param, values));
   }

   smStatus smOptimizer::setObjectiveFunction(SimTK::objectiveFunction *objFunc) {

       return(((optimizerImplementation *)data)->setObjectiveFunction(objFunc));
   }

   smStatus smOptimizer::optimize(SimTK::Vector_<Real>   &results) {
       return( ((optimizerImplementation *)data)->optimize(results));
   }

} // namespace SimTK
