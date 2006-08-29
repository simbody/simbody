
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


smHandle smMallocOptimizer( int dimension, smStatus*status){

    smHandle handle;

    SimTK::OptimizerImplementation *opt = new SimTK::OptimizerImplementation(dimension);
// TODO CATCH error from constructor
    *status = SUCCESS;
    return( (smHandle)opt);
}

smStatus  smSetCostFunction(  smHandle handle,  void (*costFunction)(int, double*,double*,double*,void*) ) {

    return( ((SimTK::OptimizerImplementation *)handle)->setObjectiveFunction(costFunction));
}

smStatus  smRunOptimizer(  smHandle handle, double *results ) {


    return( ((SimTK::OptimizerImplementation *)handle)->optimize(results));
}
    
void smFreeOptimizer(smHandle handle){
  
   delete ((SimTK::OptimizerImplementation *)handle);

   return;

}
smStatus smGetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

    return( ((SimTK::OptimizerImplementation *)handle)->getOptimizerParameters(parameter,values));
}

smStatus smSetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

  return( ((SimTK::OptimizerImplementation *)handle)->setOptimizerParameters(parameter,values));
}

smStatus smDumpOptimizerState(smHandle handle) {

     smStatus status = SUCCESS;

   return(status); 
}
