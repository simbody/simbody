/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

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
