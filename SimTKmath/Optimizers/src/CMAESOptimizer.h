#ifndef SimTK_SIMMATH_CMAES_OPTIMIZER_H_
#define SimTK_SIMMATH_CMAES_OPTIMIZER_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Michael Sherman, Nikolaus Hansen                             *
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

#include "simmath/internal/common.h"
#include "simmath/internal/OptimizerRep.h"
#include "c-cmaes/cmaes_interface.h"

#include <memory>

#if SimTK_SIMMATH_MPI
    #include <mpi.h>
    #define SimTK_COMM_WORLD MPI_COMM_WORLD
#endif

namespace SimTK {

class CMAESOptimizer: public Optimizer::OptimizerRep {
public:
    CMAESOptimizer(const OptimizerSystem& sys);
    OptimizerRep* clone() const override;
    Real optimize(SimTK::Vector& results) override;
    OptimizerAlgorithm getAlgorithm() const override { return CMAES; }

private:

    // Helper methods.
    // ---------------
    // If not using MPI, this method performs the optimization. If using MPI,
    // this delegates the evaluation of the objective function to the worker
    // processes.
    Real master(SimTK::Vector& results, const bool& useMPI);

    void checkInitialPointIsFeasible(const SimTK::Vector& x) const;

    // Wrapper around cmaes_init.
    double* init(cmaes_t& evo, Vector& results) const;
    // Edit settings in evo.sp (cmaes_readpara_t).
    void process_readpara_settings(cmaes_t& evo) const;

    void resampleToObeyLimits(cmaes_t& evo, double*const* pop);

    // May use threading or MPI.
    void evaluatePopulation(const int& popsize,
            double*const* pop, double* funvals,
            ParallelExecutor* executor, const int& nNodes);

    // ParallelExecutor Task.
    // ----------------------
    // Calls the objective function for the i-th member of the population.
    class Task : public SimTK::ParallelExecutor::Task {
    public:
        Task(CMAESOptimizer& rep, int n, double*const* pop, double* funvals)
            :   rep(rep), n(n), pop(pop), funvals(funvals) {}
        void execute(int i) override
        { rep.objectiveFuncWrapper(n, pop[i], true, &funvals[i], &rep); }
    private:
        CMAESOptimizer& rep;
        int n;
        double*const* pop;
        double* funvals;
    };

    // MPI helper methods.
    // -------------------
    #if SimTK_SIMMATH_MPI
        // If we are in a worker process (rank > 0), we just evaluate this
        // method. This method evaluates the objective function in a loop.
        Real mpi_worker(SimTK::Vector& params, const int& myRank);
        void mpi_master_sendJob(double*const* pop, const int& nParams,
                const int& ijob, const int& workerRank);
        // TODO can this be a global method I define in the cpp, so that I
        // am not polluting a header file?
        void mpi_master_receiveEval(double& buffer, MPI_Status& status,
                double* funvals);
    #endif

};

} // namespace SimTK

#endif // SimTK_SIMMATH_CMAES_OPTIMIZER_H_
