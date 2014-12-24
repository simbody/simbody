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
 * Contributors: Michael Sherman, Nikolaus Hansen, Jack Wang                  *
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

#include "CMAESOptimizer.h"

#include <bitset>

namespace SimTK {

#define SimTK_CMAES_PRINT(diag, cmds) \
do { \
if (std::bitset<2>(diag).test(0)) { printf("[SimTK CMAES] "); cmds; } \
} \
while(false)

#define SimTK_CMAES_FILE(diag, cmds) \
do { \
if (std::bitset<2>(diag).test(1)) { cmds; } \
} \
while(false)

// TODO when Simbody switches to C++11, only use unique_ptr.
#define SimTK_CMAES_SMART_PTR std::unique_ptr
#if (defined(__GNUG__) && __cplusplus<201103L)
    #undef SimTK_CMAES_SMART_PTR
    #define SimTK_CMAES_SMART_PTR std::auto_ptr
#endif

CMAESOptimizer::CMAESOptimizer(const OptimizerSystem& sys) : OptimizerRep(sys)
{
    SimTK_VALUECHECK_ALWAYS(2, sys.getNumParameters(), INT_MAX, "nParameters",
            "CMAESOptimizer");
}

Optimizer::OptimizerRep* CMAESOptimizer::clone() const {
    return new CMAESOptimizer(*this);
}

Real CMAESOptimizer::optimize(SimTK::Vector& results)
{
    // Use MPI?
    // --------
    std::string parallel;
    bool useMPI = false;
    if (getAdvancedStrOption("parallel", parallel) && (parallel == "mpi")) {
        useMPI = true;

        #if (!SimTK_SIMMATH_MPI)
            SimTK_ERRCHK_ALWAYS(false, "Optimizer::setAdvancedStrOption",
                    "You requested using MPI, but Simbody has not been "
                    "compiled with MPI. See README.md for instructions.");
        #endif
    }

    #if SimTK_SIMMATH_MPI
        if (useMPI) {
            
            // Ensure that the user has initialized MPI.
            int initialized;
            MPI_Initialized(&initialized);
            SimTK_ERRCHK_ALWAYS(initialized, "Optimizer::optimize()",
                "You requested using MPI but you did not initialize MPI. "
                "Call MPI_Init() some time before calling optimize().");

            // Are we the master node or a worker?
            int nNodes = 0;
            int myRank = 0;
            MPI_Comm_size(SimTK_COMM_WORLD, &nNodes);
            MPI_Comm_rank(SimTK_COMM_WORLD, &myRank);

            SimTK_CMAES_PRINT(diagnosticsLevel,
                    printf("Using MPI. Rank: %d.\n", myRank));

            // This is the master node
            // -----------------------
            if (myRank == 0) {

                // Run the optimization.
                Real f = master(results, useMPI);

                // Send the results to the worker nodes.
                // -------------------------------------
                // See mpi_worker for where these sends are received.
                // Send with tag 0, the end tag, to interrupt the worker loop.
                for (int workerRank = 1; workerRank < nNodes; ++workerRank) {

                    // Send the parameters.
                    MPI_Send(&results[0], results.size(), MPI_DOUBLE,
                             workerRank, /* tag: */ 0, SimTK_COMM_WORLD);
                    // Send the objective function value.
                    MPI_Send(&f, 1, MPI_DOUBLE,
                             workerRank, /* tag: */ 0, SimTK_COMM_WORLD);
                }
                return f;
            }

            // This is a worker node.
            // ----------------------
            else {
                return mpi_worker(results, myRank);
            }
        }

    #endif
    
    // Not using MPI (it's not compiled in, or it is but useMPI is false).
    // ===================================================================
    return master(results, false);
}

Real CMAESOptimizer::master(SimTK::Vector& results, const bool& useMPI)
{ 

    const OptimizerSystem& sys = getOptimizerSystem();
    int nParams = sys.getNumParameters();

    // Initialize objective function value and cmaes data structure.
    // =============================================================
    Real f; 
    cmaes_t evo;

    // Check that the initial point is feasible.
    // =========================================
    checkInitialPointIsFeasible(results);
    
    // Initialize cmaes.
    // =================
    double* funvals = init(evo, results);
    SimTK_CMAES_PRINT(diagnosticsLevel, printf("%s\n", cmaes_SayHello(&evo)));
    int popsize = cmaes_Get(&evo, "popsize");
    
    // Initialize multithreading, if requested.
    // ========================================
    std::string parallel;
    SimTK_CMAES_SMART_PTR<ParallelExecutor> exec;
    if (getAdvancedStrOption("parallel", parallel) &&
            parallel == "multithreading") {

        // Number of parallel threads.
        int nthreads = ParallelExecutor::getNumProcessors();
        getAdvancedIntOption("nthreads", nthreads);

        SimTK_CMAES_PRINT(diagnosticsLevel,
                printf("Executing on %d threads.\n", nthreads));

        exec.reset(new ParallelExecutor(nthreads));
    }

    // MPI.
    // ----
    // If using MPI, find the number of nodes / processors.
    int nNodes = 0;
    #if SimTK_SIMMATH_MPI
        if (useMPI) {
            MPI_Comm_size(SimTK_COMM_WORLD, &nNodes);
        }
    #endif
    
    // Optimize.
    // =========
    while (!cmaes_TestForTermination(&evo)) {

        // Sample a population.
        // ====================
        double*const* pop = cmaes_SamplePopulation(&evo);

        SimTK_CMAES_PRINT(diagnosticsLevel,
                printf("Generation %d.\n", (int)cmaes_Get(&evo, "generation")));

        // Resample to keep population within limits.
        // ==========================================
        resampleToObeyLimits(evo, pop);

        // Evaluate the objective function on the samples.
        // ===============================================
        evaluatePopulation(popsize, pop, funvals, exec.get(), nNodes);
        
        // Update the distribution (mean, covariance, etc.).
        // =================================================
        cmaes_UpdateDistribution(&evo, funvals);
    }

    // Wrap up.
    // ========
    SimTK_CMAES_PRINT(diagnosticsLevel,
            printf("Stop: %s\n", cmaes_TestForTermination(&evo)));

    // Update results and objective function value.
    const double* xbestever = cmaes_GetPtr(&evo, "xbestever");
    for (int i = 0; i < nParams; i++) {
        results[i] = xbestever[i]; 
    }
    f = cmaes_Get(&evo, "fbestever");

    SimTK_CMAES_FILE(diagnosticsLevel,
            cmaes_WriteToFile(&evo, "all", "allcmaes.dat"));

    // Free memory.
    cmaes_exit(&evo);
    
    return f;  
}

void CMAESOptimizer::checkInitialPointIsFeasible(const Vector& x) const {

    const OptimizerSystem& sys = getOptimizerSystem();

    if( sys.getHasLimits() ) {
        Real *lower, *upper;
        sys.getParameterLimits( &lower, &upper );
        for (int i = 0; i < sys.getNumParameters(); i++) {
            SimTK_ERRCHK4_ALWAYS(
                    lower[i] <= x[i] && x[i] <= upper[i],
                    "Optimizer::optimize",
                    "Initial guess results[%d] = %f "
                    "is not within limits [%f, %f].",
                    i, x[i], lower[i], upper[i]);
        }
    }
}

double* CMAESOptimizer::init(cmaes_t& evo, SimTK::Vector& results) const
{
    const OptimizerSystem& sys = getOptimizerSystem();
    int n = sys.getNumParameters();

    // Prepare to call cmaes_init_para.
    // ================================

    // popsize
    // -------
    int popsize = 0;
    getAdvancedIntOption("popsize", popsize);
    
    // init_stepsize
    // -------------
    double init_stepsize = 0;
    double* stddev = NULL;
    Vector init_stepsizeArray;
    if (getAdvancedRealOption("init_stepsize", init_stepsize)) {
        init_stepsizeArray.resize(n);
        for (int i = 0; i < n; i++) {
            init_stepsizeArray[i] = init_stepsize;
        }
        stddev = &init_stepsizeArray[0];
    }

    // seed
    // ----
    int seed = 0;
    if (getAdvancedIntOption("seed", seed)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(seed, "seed", "CMAESOptimizer::init");
    }

    // input parameter filename
    // ------------------------
    std::string input_parameter_filename = "none";
    SimTK_CMAES_FILE(diagnosticsLevel,
            input_parameter_filename = "writeonly";);

    // Call cmaes_init_para.
    // =====================
    // Here, we specify the subset of options that can be passed to
    // cmaes_init_para.
    cmaes_init_para(&evo,
            n,                 // dimension (0 to not set here)
            &results[0],       // xstart    (NULL to not set here)
            stddev,            // stddev    (NULL to not set here)
            seed,              // seed      (0 to not set here)
            popsize,           // lambda    (0 to not set here)
            input_parameter_filename.c_str() // input_parameter_filename
            ); 

    // Set settings that are usually read in from cmaes_initials.par.
    // ==============================================================
    process_readpara_settings(evo);

    // Once we've updated settings in cmaes_readpara_t,
    // finalize the initialization.
    return cmaes_init_final(&evo);
}

void CMAESOptimizer::process_readpara_settings(cmaes_t& evo) const
{
    // Termination criteria
    // ====================

    // stopMaxIter
    // -----------
    // maxIterations is a protected member variable of OptimizerRep.
    evo.sp.stopMaxIter = maxIterations;

    // stopTolFun
    // ----------
    // convergenceTolerance is a protected member variable of OptimizerRep.
    evo.sp.stopTolFun = convergenceTolerance;

    // stopMaxFunEvals
    // ---------------
    int stopMaxFunEvals;
    if (getAdvancedIntOption("stopMaxFunEvals", stopMaxFunEvals)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(stopMaxFunEvals, "stopMaxFunEvals",
                "Optimizer::setAdvancedIntOption");
        evo.sp.stopMaxFunEvals = stopMaxFunEvals;
    }

    // stopFitness
    // -----------
    double stopFitness;
    if (getAdvancedRealOption("stopFitness", stopFitness)) {
        evo.sp.stStopFitness.flg = 1;
        evo.sp.stStopFitness.val = stopFitness;
    }

    // stopTolFunHist
    // --------------
    double stopTolFunHist;
    if (getAdvancedRealOption("stopTolFunHist", stopTolFunHist)) {
        evo.sp.stopTolFunHist = stopTolFunHist;
    }

    // stopTolX
    // --------
    double stopTolX;
    if (getAdvancedRealOption("stopTolX", stopTolX)) {
        evo.sp.stopTolX = stopTolX;
    }

    // stopTolXFactor
    // --------------
    double stopTolUpXFactor;
    if (getAdvancedRealOption("stopTolUpXFactor", stopTolUpXFactor)) {
        evo.sp.stopTolUpXFactor = stopTolUpXFactor;
    }

    // maxtime
    // =======
    double maxtime;
    if (getAdvancedRealOption("maxTimeFractionForEigendecomposition", maxtime))
    {
        evo.sp.updateCmode.maxtime = maxtime;
    }
}

void CMAESOptimizer::resampleToObeyLimits(cmaes_t& evo, double*const* pop)
{
    const OptimizerSystem& sys = getOptimizerSystem();
    if( sys.getHasLimits() ) {

        Real *lower, *upper;
        sys.getParameterLimits( &lower, &upper );

        for (int i = 0; i < cmaes_Get(&evo, "popsize"); i++) {
            bool feasible = false; 
            while (!feasible) {
                feasible = true; 
                for (int j = 0; j < sys.getNumParameters(); j++) {
                    if (pop[i][j] < lower[j] || pop[i][j] > upper[j]) {
                        feasible = false; 
                        pop = cmaes_ReSampleSingle(&evo, i); 
                        break; 
                    }
                }
            }
        }
    }
}

void CMAESOptimizer::evaluatePopulation(
        const int& popsize, double*const* pop, double* funvals,
        ParallelExecutor* executor, const int& nNodes)
{
    const OptimizerSystem& sys = getOptimizerSystem();
    int nParams = sys.getNumParameters();

    // Execute on multiple threads.
    // ============================
    if (executor) {
        // See SimTK::CMAESOptimizer::Task; this calls the objective function.
        Task task(*this, nParams, pop, funvals);
        executor->execute(task, popsize);
    }

    // Execute across multiple processes (nodes).
    // ==========================================
    #if SimTK_SIMMATH_MPI
        else if (nNodes > 1) {
            // This only gets called by the master (root) node. Dispatch the
            // population among the worker processes, then run one evaluation
            // locally.
            // Note that popsize can be different from the number of nodes,
            // which prevents us from using a more elegant scatter/gather
            // implementation.
            int nWorkers = nNodes - 1;
            // To receive the worker's rank and tag.         
            MPI_Status status;
            // There's a job for each member of the population.
            int ijob = 0;
            // Incremented when sending a job, decremented when receiving.
            int nEvalsToReceive;
            
            // 1. Send one vector of parameters to each worker.
            // ------------------------------------------------
            // If more workers than jobs, don't give all workers a job.
            int nInitialJobs = std::min(nWorkers, popsize);
            int workerRank;
            for (ijob = 0; ijob < nInitialJobs; ++ijob)
            {
                // The worker ranks are 1 through nWorkers.
                workerRank = ijob + 1;
                mpi_master_sendJob(pop, nParams, ijob, workerRank);
            }
            nEvalsToReceive = nInitialJobs;
            
            // While waiting to receive, run an evaluation ourselves.
            objectiveFuncWrapper(nParams, pop[ijob], true, &funvals[ijob],
                    this);
            ijob++;

            // 2. As workers finish jobs, send off more jobs till all are sent.
            // ----------------------------------------------------------------
            // This temporarily stores incoming objective function evaluations.
            double buffer;
            while (ijob < popsize) {

                // Must receive a job before we can send off a new one.
                mpi_master_receiveEval(buffer, status, funvals);
                // Send the job to the worker that just finished a job.
                mpi_master_sendJob(pop, nParams, ijob, status.MPI_SOURCE);

                // Increment the count of sent jobs.
                ijob++;
            
                if (ijob < popsize) {
                    // While waiting to receive, run an evaluation ourselves.
                    objectiveFuncWrapper(nParams, pop[ijob], true,
                            &funvals[ijob], this);
                    ijob++;
                }
            }

            // 3. All jobs are sent, collect the remaining jobs as they finish.
            // ----------------------------------------------------------------
            while (nEvalsToReceive > 0) {
                mpi_master_receiveEval(buffer, status, funvals);
                nEvalsToReceive--;
            }
        }
    #endif

    // Execute normally (in serial).
    // =============================
    else {
        for (int i = 0; i < popsize; i++) {
            objectiveFuncWrapper(nParams, pop[i], true, &funvals[i], this);
        }
    }
}

#if SimTK_SIMMATH_MPI
void CMAESOptimizer::mpi_master_sendJob(double*const* pop, const int& nParams,
        const int& ijob, const int& workerRank)
{
    SimTK_CMAES_PRINT(diagnosticsLevel,
            printf("MPI: Dispatching job %d to worker %d.\n",
                ijob, workerRank));

    // The end tag is 0. To reserve 0, we use jobNumber = tagNumber - 1.
    MPI_Send(pop[ijob],           // data to send.
             nParams, MPI_DOUBLE, // size and type of data.
             workerRank,          // rank of destination.
             ijob + 1,            // tag, the job number + 1.
             SimTK_COMM_WORLD);
}

void CMAESOptimizer::mpi_master_receiveEval(double& buffer,
        MPI_Status& status, double* funvals)
{
    MPI_Recv(&buffer,        // objective value to receive.
             1, MPI_DOUBLE,  // size and type of data.
             MPI_ANY_SOURCE, // accept any worker.
             MPI_ANY_TAG,    // accept any job.
             SimTK_COMM_WORLD, &status);

    // Store the objective function value from this worker.
    // Jobs may come back in a different order than they were
    // dispatched. We use the tag (job number) to deal with this.
    int ijob = status.MPI_TAG - 1;
    funvals[ijob] = buffer;

    SimTK_CMAES_PRINT(diagnosticsLevel,
            printf("MPI: Received    job %d's evaluation. f = %f.\n",
                ijob, buffer));
}

Real CMAESOptimizer::mpi_worker(SimTK::Vector& params,
        const int& myRank) {
    const OptimizerSystem& sys = getOptimizerSystem();
    int nParams = sys.getNumParameters();

    MPI_Status status;

    // This stores objective function value.
    Real f;

    // Receive jobs, evaluate the objective function, send back the result.
    while (true) {

        // Receive a job.
        // --------------
        MPI_Recv(&params[0],          // data to send.
                 nParams, MPI_DOUBLE, // size and type of data.
                 0,                   // rank of sender (master, root).
                 MPI_ANY_TAG,         // tag.
                 SimTK_COMM_WORLD, &status);

        // Has the master node finished optimizing? Look for end tag 0.
        // ------------------------------------------------------------
        if (status.MPI_TAG == 0) {

            SimTK_CMAES_PRINT(diagnosticsLevel,
                    printf("MPI: Worker %d received end tag.\n", myRank));

            // params now contains the optimization solution/results.
            // Receive the objective function value.
            MPI_Recv(&f, 1, MPI_DOUBLE,
                     0, 0, SimTK_COMM_WORLD, MPI_STATUS_IGNORE);

            return f;
        }

        // Evaluate objective function.
        // ----------------------------
        objectiveFuncWrapper(nParams, &params[0], true, &f, this);

        // Send the objective function value back to the master.
        // -----------------------------------------------------
        MPI_Send(&f,             // objective function value to send.
                 1, MPI_DOUBLE,  // size and type of data.
                 0,              // rank of destination (master, root).
                 status.MPI_TAG, // send back with same tag as received.
                 SimTK_COMM_WORLD);
    }
}
#endif

#undef SimTK_CMAES_PRINT
#undef SimTK_CMAES_FILE
#undef SimTK_CMAES_SMART_PTR

} // namespace SimTK
