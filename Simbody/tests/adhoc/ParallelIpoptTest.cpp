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

/*This test is the parallelized version of IpoptTest in the regression suit.
* Multiple instances of Ipopt are executed concurrently on different threads
* to verfy for Ipopt thread-safety. If Ipopt is not thread-safe, there will 
* be a segmentation fault.*/

#include "SimTKmath.h"

#include <iostream>
#include <future>
#include <vector>
using std::cout;
using std::endl;
using std::future;
using std::vector;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::Real;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;


static int  NUMBER_OF_PARAMETERS = 4;
static int  NUMBER_OF_EQUALITY_CONSTRAINTS = 1;
static int  NUMBER_OF_INEQUALITY_CONSTRAINTS = 1;
static int  NUMBER_OF_THREADS = 10;
/*
* Adapted from Ipopt's hs071 example problem
*
*     min   x1*x4*(x1 + x2 + x3)  +  x3
*     s.t.  x1*x2*x3*x4                   >=  25
*           x1**2 + x2**2 + x3**2 + x4**2  =  40
*           1 <=  x1,x2,x3,x4  <= 5
*
*     Starting point:
*        x = (1, 5, 5, 1)
*
*     Optimal solution:
*        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
*
*/

class ProblemSystem : public OptimizerSystem {
public:


    int objectiveFunc(const Vector &coefficients, bool new_coefficients, Real& f) const override {
        const Real *x;

        x = &coefficients[0];

        f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
        return(0);
    }

    int gradientFunc(const Vector &coefficients, bool new_coefficients, Vector &gradient) const override {
        const Real *x;

        x = &coefficients[0];

        gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        gradient[1] = x[0] * x[3];
        gradient[2] = x[0] * x[3] + 1;
        gradient[3] = x[0] * (x[0] + x[1] + x[2]);

        return(0);

    }
    int constraintFunc(const Vector &coefficients, bool new_coefficients, Vector &constraints)  const override {
        const Real *x;

        x = &coefficients[0];
        constraints[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 40.0;
        constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;

        return(0);
    }

    int constraintJacobian(const Vector& coefficients, bool new_coefficients, Matrix& jac)  const override {
        //  int constraintJacobian( const Vector& coefficients, const bool new_coefficients, Vector& jac)  const{
        const Real *x;

        x = &coefficients[0];

        jac(0, 0) = 2 * x[0];
        jac(0, 1) = 2 * x[1];
        jac(0, 2) = 2 * x[2];
        jac(0, 3) = 2 * x[3];
        jac(1, 0) = x[1] * x[2] * x[3];
        jac(1, 1) = x[0] * x[2] * x[3];
        jac(1, 2) = x[0] * x[1] * x[3];
        jac(1, 3) = x[0] * x[1] * x[2];

        return(0);
    }

    /*   ProblemSystem() : OptimizerSystem( NUMBER_OF_PARAMETERS, NUMBER_OF_CONSTRAINTS ) {} */

    ProblemSystem(const int numParams, const int numEqualityConstraints, const int numInequalityConstraints) :
        OptimizerSystem(numParams)
    {
        setNumEqualityConstraints(numEqualityConstraints);
        setNumInequalityConstraints(numInequalityConstraints);
    }

};

struct OptimizationResult {
    Real f;
    Vector results;
    bool success;
};

void printResults(const OptimizationResult& optimizationResults) {
    printf("IpoptTest.cpp: f = %f params = ", optimizationResults.f);
    for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
        printf(" %f", optimizationResults.results[i]);
    }
    printf("\n");
}

int printError(const OptimizationResult& optimizationResults) {

    int returnValue = 0;
    static const Real TOL = 1e-4;
    Real expected[] = { 1.00000000, 4.74299963, 3.82114998, 1.37940829 };
    for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
        if (optimizationResults.results[i] > expected[i] + TOL
            || optimizationResults.results[i] < expected[i] - TOL) {
            printf(" IpoptTest.cpp: error results[%d] = %f  expected=%f \n", i,
                optimizationResults.results[i], expected[i]);
            returnValue = 1;
        }
    }
    return returnValue;
}

int main() {

    int i;
    int returnValue = 0;
    /* create the system to be optimized */
    ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_EQUALITY_CONSTRAINTS, NUMBER_OF_INEQUALITY_CONSTRAINTS);

    Vector results(NUMBER_OF_PARAMETERS);
    Vector lower_bounds(NUMBER_OF_PARAMETERS);
    Vector upper_bounds(NUMBER_OF_PARAMETERS);

    /* set initial conditions */
    results[0] = 1.0;
    results[1] = 5.0;
    results[2] = 5.0;
    results[3] = 1.0;

    /* set bounds */
    for (i = 0; i<NUMBER_OF_PARAMETERS; i++) {
        lower_bounds[i] = 1.0;
        upper_bounds[i] = 5.0;
    }

    sys.setParameterLimits(lower_bounds, upper_bounds);

    //using future to easily collect the results of the optimizations
    vector<future<OptimizationResult>> futures;

    /* compute multiple optimizations. 
    use an internal copy of opt as Optimizer cannot be 
    shared across threads.*/
    auto func([&sys](const Vector& startingValues) {
        OptimizationResult optimizationResult;
        optimizationResult.success = true;
        try {
            Optimizer opt(sys);
            opt.setConvergenceTolerance(1e-4);
            opt.setDiagnosticsLevel(0);
            opt.setLimitedMemoryHistory(500); // works well for our small systems
            opt.setAdvancedBoolOption("warm_start", true);
            opt.setAdvancedRealOption("obj_scaling_factor", 1);
            opt.setAdvancedRealOption("nlp_scaling_max_gradient", 1);
            optimizationResult.results = startingValues;
            optimizationResult.f = opt.optimize(optimizationResult.results);
        }
        catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            optimizationResult.success = false; // failure
            printf("IpoptTest.cpp: Caught exception \n");
        }
        return optimizationResult;
    });

    //launch threads
    for (int iThread(0); iThread < NUMBER_OF_THREADS; ++iThread) {
        futures.emplace_back(
            std::async(std::launch::async, func, results));
    }
 
    //collect optimization results
    cout << "Printing the result of " << NUMBER_OF_THREADS 
        << " concurrent optimizations.\n";
    for (auto& f : futures) {
        f.wait();
        auto ans(f.get());
        returnValue |= (!ans.success);
        printResults(ans);
        returnValue |= printError(ans);
    }

    return(returnValue);
}
