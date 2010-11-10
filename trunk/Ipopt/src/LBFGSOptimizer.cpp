

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
#include "LBFGSOptimizer.h"

using std::cout;
using std::endl;

namespace SimTK {

Optimizer::OptimizerRep* LBFGSOptimizer::clone() const {
    return( new LBFGSOptimizer(*this) );
}


LBFGSOptimizer::LBFGSOptimizer( const OptimizerSystem& sys )
    : OptimizerRep( sys ) ,
      xtol(1e-16)

{
     /* internal flags for LBFGS */
     iprint[0] = iprint[1] = 0; // no output generated

     if( sys.getNumParameters() < 1 ) {
        const char* where = "Optimizer Initialization";
        const char* szName = "dimension";
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  sys.getNumParameters(), INT_MAX, where); 
     }
} 

Real LBFGSOptimizer::optimize(  Vector &results ) {
    int iflag[1] = {0};
    Real f;
    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();
    int m = limitedMemoryHistory;

    iprint[0] = iprint[1] = iprint[2] = diagnosticsLevel; 

    double tol;
    if( getAdvancedRealOption("xtol", tol ) ) {
        SimTK_APIARGCHECK_ALWAYS(tol > 0,"LBFGSOptimizer","optimize",
        "xtol must be positive \n");
        xtol = tol;
    }

    lbfgs_( n, m, &results[0], &f, iprint, &convergenceTolerance,  &xtol );

    objectiveFuncWrapper(n, &results[0], true, &f, this);
    return f;
}

} // namespace SimTK
