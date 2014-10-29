/* -------------------------------------------------------------------------- *
 *     Simbody(tm) Example: Global optimization with CMAES                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

#include "SimTKmath.h"
#include <iostream>

// If you compiled Simbody with CMake variable SIMBODY_ENABLE_MPI set to ON,
// try changing the value of USE_MPI to 1. If you don't, then the optimization
// will run completely separately on each MPI process.
#define USE_MPI 0

// Run this example with a command like $ mpirun -n 8 ./CMAESOptimization

#if USE_MPI
    #include <mpi.h>
#endif

using namespace SimTK;

/* A very complex 2D function with many local minima.
 * http://www.sfu.ca/~ssurjano/drop.html
 *
 *     min   -(1 + cos(12 * (x1^2 + x2^2)^(1/2))) / (0.5 * (x1^2 + x2^2) + 2)
 *     s.t.  -5.12 <= x1,x2 <= 5.12
 *
 *     Starting point:
 *        x = (2, 2)
 *
 *     Optimal solution:
 *        x = (0, 0)
 *
 */
class DropWave : public OptimizerSystem {
public:
    DropWave() : OptimizerSystem(2) {
        // The website above says this function usually has the following
        // bounds:
        Vector limits(getNumParameters());
        limits.setTo(5.12);
        setParameterLimits(-limits, limits);
    }
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        const Real dotprod = x[0] * x[0] + x[1] * x[1];
        f = -(1 + cos(12 * sqrt(dotprod))) / (0.5 * dotprod + 2);
        return 0;
    }
};


int main(int argc, char* argv[]) {

    #if USE_MPI
        MPI_Init(&argc, &argv);
    #endif

    // Create the system to be optimized.
    DropWave sys;

    int nParameters = sys.getNumParameters();

    // Set initial guess.
    Vector results(nParameters);
    results.setTo(2);

    // Create the Optimizer.
    Optimizer opt(sys, SimTK::CMAES); 

    // Configure the optimizer.
    opt.setConvergenceTolerance(1e-5);
    // opt.setDiagnosticsLevel(3);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("sigma", 3.5);
    opt.setAdvancedIntOption("lambda", 1000);
    opt.setAdvancedIntOption("stopMaxFunEvals", 100000);
    opt.setAdvancedIntOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    // opt.setAdvancedStrOption("parallel", "mulithreading");

    #if USE_MPI
        opt.setAdvancedStrOption("parallel", "mpi");
    #endif

    Real f = opt.optimize(results);

    #if USE_MPI
        int myRank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        // We only want the master / root node to print the result.
        if (myRank == 0) {
    #endif
    printf("Optimal Solution: f = %f   parameters = %f %f \n",
            f, results[0], results[1]);
    #if USE_MPI
        }
    #endif  

    #if USE_MPI
        MPI_Finalize();
    #endif

    return 0;
}
