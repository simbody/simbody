
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


/* 
???? separtate classes for each type of optimizer or one class that chooses the
     appropriate algorithm ?

*/

// #include "objectiveFunction.h"

#include "optimizerInterface.h"
#include "optimizerImplementation.h"

namespace SimTK {

// using namespace SimTK; TODO resolve this

class smOptimizer : public smObject,  public smOptimizerInterface {
   public:
    smOptimizer( int );

int  setOptimizerParameters(int param, double *values) {

         return(((optimizerImplementation *)data)->setOptimizerParameters(param, values));
      }

int getOptimizerParameters(int param, double *values) {

         return(((optimizerImplementation *)data)->getOptimizerParameters(param, values));
     }

int setObjectiveFunction(int *func) {

       return(((optimizerImplementation *)data)->setObjectiveFunction(func));
}
int setInitialPosition(double *pos) {

       return(((optimizerImplementation *)data)->setInitialPosition(pos));
}
int optimize(double *results) {

       return( ((optimizerImplementation *)data)->optimize(results));
}

}; // Class smOptimizer
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_H

