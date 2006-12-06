#ifndef _SimTK_SIMMATH_H_
#define _SimTK_SIMMATH_H_

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

typedef void* smHandle; 

enum { 	OPTIMIZER, ROOT_FINDER };
typedef enum { LBFGS, LBFGSB,  InteriorPoint } OptimizerAlgorithm;
enum { SUCCESS, MALLOC_FAILED, UNKNOWN_OPTIMIZER, UNKNOWN_PARAMETER, INVALID_VALUE };
enum { TRACE, MAX_FUNCTION_EVALUATIONS, DEFAULT_STEP_LENGTH, LINE_SEARCH_ACCURACY, GRADIENT_CONVERGENCE_TOLERANCE };


#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "Simmath/Optimizer.h"

namespace SimTK {

namespace Exception {

class OptimizerFailed : public Base {
public:
        OptimizerFailed( const char * fn, int ln, String msg) : Base(fn, ln)
        {
            setMessage("Optmizer failed: " + msg );
        }
private:
};

class UnrecognizedParameter : public Base {
public:
        UnrecognizedParameter( const char * fn, int ln, String msg) : Base(fn, ln)
        {
            setMessage("Unrecognized Parameter: " + msg );
        }
private:
};

} // namespace Exception

} //  namespace SimTK


#endif //_SimTK_SIMMATH_H_
