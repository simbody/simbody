#ifndef _SimTK_OPTIMIZER_INTERFACE_H_
#define _SimTK_OPTIMIZER_INTERFACE_H_

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
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "ObjectiveFunction.h"


namespace SimTK {

class OptimizerInterface { 

public:
    virtual ~OptimizerInterface() {};
   

    /* by default class chooses algorithm based on size of problem. 
       The algorithm can be overridden by setting the Algo param. */
       /* should there be a set method for each parameter ? */
    virtual void setOptimizerParameters(unsigned int, double *) = 0;
    virtual void getOptimizerParameters(unsigned int, double *) = 0;
/*
    virtual void getOptimizerParameters(unsigned int, int *) = 0;
    virtual void setOptimizerParameters(unsigned int, int *) = 0;
*/

    virtual double optimize(SimTK::Vector &) = 0; // checks to see if space needs to be allocated

       // constructor sets internal flag to require allocation
       // if dimension or algorithm changes the reallocate flag is
       // set to trigger freeing the old mem and allocating new.


}; // end class optimizeInterface
} // namespace SimTK
#endif  //_SimTK_OPTIMIZER_INTERFACE_H_
