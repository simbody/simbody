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

#include "InteriorPointOptimizer.h"

using std::cout;
using std::endl;


namespace SimTK {

static std::string applicationReturnStatusToString(int status) {
    switch(status) {
        case Solve_Succeeded: return "Solve succeeded";
        case Solved_To_Acceptable_Level: return "Solved to acceptable level";
        case Infeasible_Problem_Detected: return "Infeasible problem detected";
        case Search_Direction_Becomes_Too_Small: return "Search direction becomes too small";
        case Diverging_Iterates: return "Diverging iterates";
        case User_Requested_Stop: return "User requested stop";
        case Maximum_Iterations_Exceeded: return "Maximum iterations exceeded";
        case Restoration_Failed: return "Restoration failed";
        case Error_In_Step_Computation: return "Error in step computation";
        case Not_Enough_Degrees_Of_Freedom: return "Not enough degrees of freedom";
        case Invalid_Problem_Definition: return "Invalid problem definition";
        case Invalid_Option: return "Invalid option";
        case Invalid_Number_Detected: return "Invalid number detected";
        case Unrecoverable_Exception: return "Unrecoverable exception";
        case NonIpopt_Exception_Thrown: return "Non-Ipopt exception thrown";
        case Insufficient_Memory: return "Insufficient memory";
        case Internal_Error: return "Internal error";
        default: return "Unknown Ipopt return status";
    }
}
Optimizer::OptimizerRep* InteriorPointOptimizer::clone() const {
	return( new InteriorPointOptimizer(*this) );
}

// Assume by the time this constructor is called, the number of parameters and constraints has been finalized
InteriorPointOptimizer::InteriorPointOptimizer( const OptimizerSystem& sys )
        : OptimizerRep( sys ) {

        int n = sys.getNumParameters();
        int m = sys.getNumConstraints();

        if( n < 1 ) {
            const char* where = " InteriorPointOptimizer Initialization";
            const char* szName = "dimension";
            SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1, n, INT_MAX, where); 
        }

        // Initialize arrays to store multipliers -- will be used for warm starts
        mult_x_L = new Number[n];
        mult_x_U = new Number[n];
        mult_g = new Number[m];
        for(int i=0;i<n;i++) mult_x_L[i] = mult_x_U[i] = 0;
        for(int i=0;i<m;i++) mult_g[i] = 0;

        g_L = new Number[m];
        g_U = new Number[m];
        /* set the bounds on the equality constraint functions */
        for(int i=0;i<sys.getNumEqualityConstraints();i++){
            g_U[i] = g_L[i] = 0.0;
        }
        /* set the bounds on the inequality constraint functions */
        for(int i=sys.getNumEqualityConstraints();i<m;i++){
            g_U[i] = SimTK::Real(POSITIVE_INF);
            g_L[i] = 0.0;
        }

        firstOptimization = true;
    } 


    SimTK::Real InteriorPointOptimizer::optimize(  Vector &results ) {

        int n = getOptimizerSystem().getNumParameters();
        int m = getOptimizerSystem().getNumConstraints();

        Index index_style = 0; /* C-style; start counting of rows and column indices at 0 */
        Index nele_hess = 0;
        Index nele_jac = n*m; /* always assume dense */

        // Parameter limits
        Number *x_L = NULL, *x_U = NULL;
        if( getOptimizerSystem().getHasLimits() ) {
           getOptimizerSystem().getParameterLimits( &x_L, &x_U);
        } else {
           x_U = new Number[n];
           x_L = new Number[n];
           for(int i=0;i<n;i++) {
              x_U[i] = SimTK::Real(POSITIVE_INF);
              x_L[i] = SimTK::Real(NEGATIVE_INF);
           }
        }

        SimTK::Real *x = &results[0];

        IpoptProblem nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 
                           nele_hess, index_style, objectiveFuncWrapper, constraintFuncWrapper, 
                           gradientFuncWrapper, constraintJacobianWrapper, hessianWrapper);

        // If you want to verify which options are getting set in the optimizer, you can create a file ipopt.opt
        // with "print_user_options yes", and set print_level to (at least 1).  It will then print the options to the screen.
        
        // sherm 100302: you have to set all of these tolerances to get IpOpt to change
        // its convergence criteria; see OptimalityErrorConvergenceCheck::CheckConvergence().
        // We'll set acceptable tolerances to the same value to disable them.
        AddIpoptNumOption(nlp, "tol", convergenceTolerance);
        AddIpoptNumOption(nlp, "dual_inf_tol", convergenceTolerance);
        AddIpoptNumOption(nlp, "constr_viol_tol", constraintTolerance);
        AddIpoptNumOption(nlp, "compl_inf_tol", convergenceTolerance);
        AddIpoptNumOption(nlp, "acceptable_tol", convergenceTolerance);
        AddIpoptNumOption(nlp, "acceptable_dual_inf_tol", convergenceTolerance);
        AddIpoptNumOption(nlp, "acceptable_constr_viol_tol", constraintTolerance);
        AddIpoptNumOption(nlp, "acceptable_compl_inf_tol", convergenceTolerance);


        AddIpoptIntOption(nlp, "max_iter", maxIterations);
        AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
        AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory"); // needs to be limited-memory unless you have explicit hessians
        AddIpoptIntOption(nlp, "limited_memory_max_history", limitedMemoryHistory);
        AddIpoptIntOption(nlp, "print_level", diagnosticsLevel); // default is 4

        int i;
        static const char *advancedRealOptions[] = {
                                                    "compl_inf_tol", 
                                                    "dual_inf_tol", 
                                                    "constr_viol_tol", 
                                                    "acceptable_tol", 
                                                    "acceptable_compl_inf_tol", 
                                                    "acceptable_constr_viol_tol", 
                                                    "acceptable_dual_inf_tol", 
                                                    "diverging_iterates_tol", 
                                                    "barrier_tol_factor", 
                                                    "obj_scaling_factor", 
                                                    "nlp_scaling_max_gradient", 
                                                    "bounds_relax_factor", 
                                                    "recalc_y_feas_tol", 
                                                    "mu_init", 
                                                    "mu_max_fact", 
                                                    "mu_max", 
                                                    "mu_min", 
                                                    "mu_linear_decrease_factor", 
                                                    "mu_superlinear_decrease_factor", 
                                                    "bound_frac", 
                                                    "bound_push", 
                                                    "bound_mult_init_val", 
                                                    "constr_mult_init_max", 
                                                    "constr_mult_init_val", 
                                                    "warm_start_bound_push", 
                                                    "warm_start_bound_frac", 
                                                    "warm_start_mult_bound_push", 
                                                    "warm_start_mult_init_max", 
                                                    "recalc_y_feas_tol", 
                                                    "expect_infeasible_problem_ctol", 
                                                    "soft_resto_pderror_reduction_factor", 
                                                    "required_infeasibility_reduction", 
                                                    "bound_mult_reset_threshold", 
                                                    "constr_mult_reset_threshold", 
                                                    "max_hessian_perturbation", 
                                                    "min_hessian_perturbation", 
                                                    "first_hessian_perturbation", 
                                                    "perturb_inc_fact_first", 
                                                    "perturb_inc_fact", 
                                                    "perturb_dec_fact", 
                                                    "jacobian_reqularization_value", 
                                                    "derivative_test_perturbation", 
                                                    "derivative_test_tol", 
                                                    0}; 
        Real value;
        for(i=0;advancedRealOptions[i];i++) {
            if(getAdvancedRealOption(advancedRealOptions[i],value))
                AddIpoptNumOption(nlp, advancedRealOptions[i], value);
        }

        static const std::string advancedStrOptions[] = {"nlp_scaling_method",
                                                         "honor_original_bounds", 
                                                         "check_derivatives_for_naninf", 
                                                         "mu_strategy", 
                                                         "mu_oracle", 
                                                         "fixed_mu_oracle", 
                                                         "alpha_for_y", 
                                                         "recalc_y", 
                                                         "expect_infeasible_problem", 
                                                         "print_options_documentation", 
                                                         "print_user_options", 
                                                         "start_with_resto", 
                                                         "evaluate_orig_obj_at_resto_trial", 
                                                         "hessian_approximation", 
                                                         "derivative_test", 
                                                         ""}; 
        std::string svalue;
        for(i=0;!advancedStrOptions[i].empty();i++) {
            if(getAdvancedStrOption(advancedStrOptions[i], svalue))
                AddIpoptStrOption(nlp, advancedStrOptions[i].c_str(), svalue.c_str());
        }

        static const char*  advancedIntOptions[] = {"quality_function_max_section_steps",
                                                         "max_soc",
                                                         "watchdog_shorted_iter_trigger",
                                                         "watchdog_trial_iter_max",
                                                         "max_refinement_steps",
                                                         "min_refinement_steps",
                                                         "limited_memory_max_history",
                                                         "limited_memory_max_skipping",
                                                         "derivative_test_print_all",
                                                         0}; 
        int ivalue;
        for(i=0;advancedIntOptions[i];i++) {
            if(getAdvancedIntOption(advancedIntOptions[i], ivalue))
                AddIpoptIntOption(nlp, advancedIntOptions[i], ivalue);
        }

        // Only makes sense to do a warm start if this is not the first call to optimize() (since we need 
        // reasonable starting multiplier values)
        bool use_warm_start=false;
        if(getAdvancedBoolOption("warm_start", use_warm_start) && use_warm_start && !firstOptimization) {
            AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
            AddIpoptStrOption(nlp, "warm_start_entire_iterate", "yes");
            //AddIpoptStrOption(nlp, "warm_start_same_structure", "yes"); // couldn't get this one to work
        } 

        SimTK::Real obj;

        int status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, (void *)this );

        FreeIpoptProblem(nlp); 

        // Only delete these if they aren't pointing to existing parameter limits
        if( !getOptimizerSystem().getHasLimits() ) {
           delete [] x_U;
           delete [] x_L;
        }

        if(status == Solved_To_Acceptable_Level) {
            std::cout << "Ipopt: Solved to acceptable level" << std::endl;
        } else if (status != Solve_Succeeded) {
            if( status != NonIpopt_Exception_Thrown) {
                char buf[1024];
                sprintf(buf, "Ipopt: %s (status %d)",applicationReturnStatusToString(status).c_str(),status);
                SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf));
            }
        }

        firstOptimization = false;

        return(obj);
    }


} // namespace SimTK
