/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
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

namespace SimTK {

CMAESOptimizer::CMAESOptimizer(const OptimizerSystem& sys) : OptimizerRep(sys),
    m_seed(0)
{
    SimTK_VALUECHECK_ALWAYS(2, sys.getNumParameters(), INT_MAX, "nParameters",
            "CMAESOptimizer");
}

Optimizer::OptimizerRep* CMAESOptimizer::clone() const {
    return( new CMAESOptimizer(*this) );
}

Real CMAESOptimizer::optimize(SimTK::Vector& results)
{ 
    Real f = HUGE_VAL; 
    const OptimizerSystem& sys = getOptimizerSystem();
    cmaes_t evo;
    int n = sys.getNumParameters();

    // TODO make sure initial point is feasible.

    // Prepare to call cmaes_init.
    // ===========================

    // TODO move to processBefore
	int numsamples = 0;
    getAdvancedIntOption("lambda", numsamples );
	if (numsamples == 0) {
		numsamples = 4+int(3*std::log((double)n)); 
        // TODO unnecessary.
	}
	
    // TODO move to processBefore, and change to 0.3
	double stepsize = 0;
    getAdvancedRealOption("sigma", stepsize ); 
	if (stepsize == 0.0) {
		stepsize = 0.1; 
	}
	Vector stepsizeArray(n);
	for (int i = 0; i < n; i++) {
		stepsizeArray[i] = stepsize;  
	}
	
    // TODO clean up.
	bool isresume = false; 
    getAdvancedBoolOption("resume", isresume );

    processSettingsBeforeCMAESInit(evo);

    // cmaes_init.
    // ===========
	 
	double* funvals = cmaes_init(&evo,
            n,                 // dimension
            &results[0],       // xstart
            &stepsizeArray[0], // stddev
            m_seed,            // seed
            numsamples,        // lambda
            "non"              // input_parameter_filename
            ); 
	if (isresume) {
		cmaes_resume_distribution(&evo, (char*)"resumecmaes.dat"); 
	}

    // TODO when to call this?
    processSettingsAfterCMAESInit(evo);

    if (diagnosticsLevel == 1) {
        printf("%s\n", cmaes_SayHello(&evo));
    }

    // TODO let the master also run an objective function evaluation.
	
    // Optimize.
    // =========
	while (!cmaes_TestForTermination(&evo)) {

        // Sample a population.
        // ====================
		double*const* pop = cmaes_SamplePopulation(&evo);

        // Resample to keep population within limits.
        // ==========================================
    	if( sys.getHasLimits() ) {
    		Real *lowerLimits, *upperLimits;
        	sys.getParameterLimits( &lowerLimits, &upperLimits );
			
			for (int i = 0; i < cmaes_Get(&evo, "popsize"); i++) {
				bool feasible = false; 
				while (!feasible) {
					feasible = true; 
					for (int j = 0; j < n; j++) {
						if (pop[i][j] < lowerLimits[j] || 
								pop[i][j] > upperLimits[j]) {
							feasible = false; 
							pop = cmaes_ReSampleSingle(&evo, i); 
							break; 
						}
					}
				}
			}
		}

        // Evaluate the objective function on the samples.
        // ===============================================
        for (int i = 0; i < cmaes_Get(&evo, "lambda"); i++) {
            objectiveFuncWrapper(n, pop[i], true, &funvals[i], this);
        }
		
        // Update the distribution (mean, covariance, etc.).
        // =================================================
		cmaes_UpdateDistribution(&evo, funvals);
		
        // Update best-yet parameters and objective function value.
        // ========================================================
		const double* optx = cmaes_GetPtr(&evo, "xbestever");
		for (int i = 0; i < n; i++) {
			results[i] = optx[i]; 
		}
		f = cmaes_Get(&evo, "fbestever");
	}

    // Wrap up.
    // ========
    printf("Stop:\n%s\n", cmaes_TestForTermination(&evo));
	cmaes_WriteToFile(&evo, "resume", "resumecmaes.dat");
    cmaes_exit(&evo);
	
	return f;  
}

void CMAESOptimizer::processSettingsBeforeCMAESInit(cmaes_t& evo)
{
    // seed
    // ====
    // If the user has provided a seed, check it.
    if (getAdvancedIntOption("seed", m_seed)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(m_seed, "seed",
                "CMAESOptimizer::processSettingsBeforeCMAESInit");
        // cmaes' seed is not set here. Seed is used in the call to cmaes_init.
    }
}

void CMAESOptimizer::processSettingsAfterCMAESInit(cmaes_t& evo)
{
    // stopMaxIter
    // ===========
    // maxIterations is a protected member variable of OptimizerRep.
    // TODO this should be set earlier, because parameters depend on what
    // stopMaxIter is.
	evo.sp.stopMaxIter = maxIterations;

    // stopTolFun
    // ==========
    evo.sp.stopTolFun = convergenceTolerance;

    // stopMaxFunEvals
    // ===============
    int stopMaxFunEvals;
    if (getAdvancedIntOption("stopMaxFunEvals", stopMaxFunEvals)) {
        SimTK_VALUECHECK_NONNEG_ALWAYS(stopMaxFunEvals, "stopMaxFunEvals",
                "CMAESOptimizer::processSettingsAfterCMAESInit");
        evo.sp.stopMaxFunEvals = stopMaxFunEvals;
    }

    // maxtime
    // =======
    double maxtime;
    if (getAdvancedRealOption("maxTimeFractionForEigendecomposition", maxtime))
    {
        evo.sp.updateCmode.maxtime = maxtime;
    }

    // TODO should call readpara_SupplementDefaults(), as this updates certain
    // parameters based on the provided values (e.g., damps).
}

} // namespace SimTK
