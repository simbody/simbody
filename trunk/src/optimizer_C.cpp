
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

namespace SimTK {

smHandle smMallocOptimizer( int dimension, int nConstraints, int nEqualConstraints, int nBounds){

    smHandle handle;

    Optimizer* opt = new Optimizer(dimension, nConstraints, nEqualConstraints, nBounds);
    return( (smHandle)opt);
}

void  smSetCostFunction(  smHandle handle,  void (*costFunction)(int, double*,double*,double*,void*) ) {

    ((SimTK::OptimizerInterface *)handle)->setObjectiveFunction(costFunction);
    return;
}

void  smRunOptimizer(  smHandle handle, double *results ) {


    ((SimTK::OptimizerInterface *)handle)->optimize(results);
    return;
}
    
void smFreeOptimizer(smHandle handle){
 
   // TODO need to call the destructor for the correct implementation 
   delete ((SimTK::OptimizerInterface *)handle);

   return;

}
void smGetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

    ((SimTK::OptimizerInterface *)handle)->getOptimizerParameters(parameter,values);
    return;
}

void smSetOptimizerParameters( smHandle handle, unsigned int parameter, double *values){

  ((SimTK::OptimizerInterface *)handle)->setOptimizerParameters(parameter,values);
  return;
}

void smDumpOptimizerState(smHandle handle) {

   return; 
}
 
}
