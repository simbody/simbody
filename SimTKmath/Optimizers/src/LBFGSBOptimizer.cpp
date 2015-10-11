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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "LBFGSBOptimizer.h"
#include <cstring>

using std::cout;
using std::endl;

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strcat, sprintf, etc.
#endif

namespace SimTK {

Optimizer::OptimizerRep* LBFGSBOptimizer::clone() const {
    return new LBFGSBOptimizer(*this);
}

LBFGSBOptimizer::LBFGSBOptimizer( const OptimizerSystem& sys )
:   OptimizerRep( sys ),
    factr( 1.0e7) {
    int n,i;

    n = sys.getNumParameters();

    if( n < 1 ) {
        const char* where = "Optimizer Initialization";
        const char* szName = "dimension";
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  n, INT_MAX, where);
    }

    /* We don't yet know what kinds of bounds we'll have so set the bounds
       descriptor nbd to an illegal value. */
    nbd = new int[n];
    for(i=0;i<n;i++)
        nbd[i] = -1;
}

Real LBFGSBOptimizer::optimize(  Vector &results ) {
    int run_optimizer = 1;
    char task[61];
    Real f;
    int *iwa;
    char csave[61];
    bool lsave[4];
    int isave[44];
    Real dsave[29];
    Real *wa;
    Real *lowerLimits, *upperLimits;
    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();
    int m = limitedMemoryHistory;
    Real *gradient;
    gradient = new Real[n];

    iprint[0] = iprint[1] = iprint[2] = diagnosticsLevel;

    if( sys.getHasLimits() ) {
        sys.getParameterLimits( &lowerLimits, &upperLimits );
        // Determine what kind of limits each parameter has.
        // nbd = 0, unbounded; 1, only lower; 2, both; 3 only upper
        for (int i=0; i < n; ++i) {
            if (lowerLimits[i] == -Infinity)
                 nbd[i] = (upperLimits[i] == Infinity ? 0 : 3);
            else nbd[i] = (upperLimits[i] == Infinity ? 1 : 2);
        }
    } else {
        lowerLimits = 0;
        upperLimits = 0;
        for(int i=0;i<n;i++)
            nbd[i] = 0;          // unbounded
    }

    iwa = new int[3*n];
    wa = new Real[((2*m + 4)*n + 12*m*m + 12*m)];

    Real factor;
    if( getAdvancedRealOption("factr", factor ) ) {
        SimTK_APIARGCHECK_ALWAYS(factor > 0,"LBFGSBOptimizer","optimize",
                                 "factr must be positive \n");
        factr = factor;
    }
    strcpy( task, "START" );
    while( run_optimizer ) {
        setulb_(&n, &m, &results[0], lowerLimits,
                upperLimits, nbd, &f, gradient,
                &factr, &convergenceTolerance, wa, iwa,
                task, iprint, csave, lsave, isave, dsave, 60, 60);

        if( strncmp( task, "FG", 2) == 0 ) {
            objectiveFuncWrapper( n, &results[0],  true, &f, this);
            gradientFuncWrapper( n,  &results[0],  false, gradient, this);
        } else if( strncmp( task, "NEW_X", 5) == 0 ){
            //objectiveFuncWrapper( n, &results[0],  true, &f, (void*)this );
        } else {
            run_optimizer = 0;
            if( strncmp( task, "CONV", 4) != 0 ){
                delete[] gradient;
                delete[] iwa;
                delete[] wa;
                SimTK_THROW1(SimTK::Exception::OptimizerFailed , SimTK::String(task) );
            }
        }
    }
    delete[] gradient;
    delete[] iwa;
    delete[] wa;

    return f;
}


} // namespace SimTK
