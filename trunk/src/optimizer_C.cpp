
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
#include "OptimizerCommon.h"
#include "Optimizer.h"

namespace SimTK {
extern "C" {
void *SimTK_mallocOptimizer( int dimension, int nConstraints, int nEqualConstraints, int nBounds){

    return( (void *)new Optimizer(dimension, nConstraints, nEqualConstraints, nBounds));
}

void  SimTK_setObjectiveFunction(  void *optimizer,  double (*f)(int, int, double*,void*) ) {

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setObjectiveFunction( f );

    return;
}

void  SimTK_setGradientFunction(  void *optimizer,  void (*f)(int, int, double*,double*,void*) ) {

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setGradientFunction( f );

    return;
}
void  SimTK_setConstraintsFunction(  void *optimizer,  void (*f)(int, int, int, double*,double*,void*) ) {

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setConstraintsFunction( f );

    return;
}
void  SimTK_setConstraintsJacobian(  void *optimizer ,  void (*f)(int, int, int, double*,double*,void*) ) {

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setConstraintsJacobian( f );

    return;
}
void SimTK_setObjectiveAndGradient( void *optimizer,  double (*f)(int, int, double*, double*, void*)){
    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setObjectiveAndGradient( f );
    return;
} 

void SimTK_setComputeHessian(  void *optimizer, void (*f)(int, int, int, double*, double*,void*) ){

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setComputeHessian( f );
    return;
} 

void  SimtK_runOptimizer(  void *optimizer, double *results ) {

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->optimize( results );
    return;
}
    
void SimTK_freeOptimizer(void *optimizer){
 
   delete ((OptimizerCommon*)optimizer);
   return;

}
void SimTK_setOptimizerParameters( void *optimizer, unsigned int parameter, double *values){

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->setOptimizerParameters(parameter,values);
    return;
}
void SimTK_getOptimizerParameters( void *optimizer, unsigned int parameter, double *values){

    (((OptimizerCommon *)((Optimizer *)optimizer)->data))->getOptimizerParameters(parameter,values);
    return;
}

} // extern "C" 
} // namespace SimTK
