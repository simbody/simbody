#include "Simmath.h"
#include "Simmath_Ctypes.h"
#include "optimizer.h"
#include "optimizerImplementation.h"


smHandle smMallocOptimizer( int dimension, smStatus*status){

    smHandle handle;

    SimTK::optimizerImplementation *opt = new SimTK::optimizerImplementation(dimension);
// TODO CATCH error from constructor
    *status = SUCCESS;
    return( (smHandle)opt);
}

smStatus  smSetCostFunction(  smHandle handle,  void (*costFunction)(double*,double*,double*) ) {

    return( ((SimTK::optimizerImplementation *)handle)->setObjectiveFunction(costFunction));
}

smStatus  smRunOptimizer(  smHandle handle, double *results ) {


    return( ((SimTK::optimizerImplementation *)handle)->optimize(results));
}
    
void smFreeOptimizer(smHandle handle){
  
   delete ((SimTK::optimizerImplementation *)handle);

   return;

}
smStatus smGetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

    return( ((SimTK::optimizerImplementation *)handle)->getOptimizerParameters(parameter,values));
}

smStatus smSetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

  return( ((SimTK::optimizerImplementation *)handle)->setOptimizerParameters(parameter,values));
}

smStatus smDumpOptimizerState(smHandle handle) {

     smStatus status = SUCCESS;

   return(status); 
}
