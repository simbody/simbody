
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
#include "OptimizerRep.h"

namespace SimTK {
    void OptimizerRep::useNumericalGradient( const bool flag ) {

       if( flag ) {     // turn on numerical gradients
           if( !numericalGradient ) { // 
              initNumericalGrad();       
              numericalGradient = true;
           }
       } else {        // turn off numerical graidents
           if( numericalGradient ) {  
              disableNumericalGrad();       
              numericalGradient = false;
           }
       }
       return;
    }
    void OptimizerRep::useNumericalJacobian( const bool flag ) {

       if( flag ) {     // turn on numerical gradients
           if( !numericalJacobian ) { // 
              initNumericalJac();       
              numericalJacobian = true;
           }
       } else {
           if( numericalJacobian ) { // 
              disableNumericalJac();       
              numericalJacobian = false;
           }
       }
       return;
    }

    void OptimizerRep::disableNumericalGrad() {  // instaniates a gradient Differentiator
       // TODO finish this 
       return;
    }
    void OptimizerRep::disableNumericalJac() {  // instaniates a jacobian Differentiator
       // TODO finish this 
       return;
    }
    void OptimizerRep::initNumericalJac() {  // instaniates a jacobian Differentiator

        cf      = new SysConstraintFunc(sysp->numParameters, sysp->numConstraints, sysp );
        jacDiff = new Differentiator(*cf);  // construct Differentiator

    }
    void OptimizerRep::initNumericalGrad() {  // instaniates a gradient Differentiator

        of       = new SysObjectiveFunc( sysp->numParameters, sysp );
        gradDiff = new Differentiator(*of );  // construct Differentiator

     }


} // namespace SimTK
