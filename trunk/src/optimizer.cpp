#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

 namespace SimTK {
smOptimizer::smOptimizer(int dimension) {
        optimizerImplementation *optPtr =  new optimizerImplementation(dimension);
        optPtr->dimension = dimension;
        data = (void *)optPtr;
}

int  smOptimizer::setOptimizerParameters(unsigned int param, double *values) {

         return(((optimizerImplementation *)data)->setOptimizerParameters(param, values));
      }
int smOptimizer::getOptimizerParameters(unsigned int param, double *values) {

         return(((optimizerImplementation *)data)->getOptimizerParameters(param, values));
     }
int smOptimizer::setObjectiveFunction(void (*func)(double*,double*,double*)) {

       return(((optimizerImplementation *)data)->setObjectiveFunction(func));
}

template < typename T >
    smStatus smOptimizer::optimize(T results) { // checks to see if space needs to be allocated
       return( ((optimizerImplementation *)data)->optimize(results));
}

template < int N, typename T, int S >
    smStatus smOptimizer::optimize(SimTK::Vec< N, T, S>  results) {
       return( ((optimizerImplementation *)data)->optimize(results));
}
template smStatus smOptimizer::optimize(double *);

template smStatus smOptimizer::optimize<2, double, 1>(Vec<2, double, 1>);
template smStatus smOptimizer::optimize(Vec<3, double, 1>);


} // namespace SimTK
