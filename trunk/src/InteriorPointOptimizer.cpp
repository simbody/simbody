

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
#include "simmath/internal/InteriorPointOptimizer.h"

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

// Assume by the time this constructor is called, the number of parameters and constraints has been finalized
InteriorPointOptimizer::InteriorPointOptimizer( OptimizerSystem& sys )
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
            g_U[i] = POSITIVE_INF;
            g_L[i] = 0.0;
        }

        firstOptimization = true;
    } 


    double InteriorPointOptimizer::optimize(  Vector &results ) {

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
              x_U[i] = POSITIVE_INF;
              x_L[i] = NEGATIVE_INF;
           }
        }

        double *x = &results[0];

        IpoptProblem nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 
                           nele_hess, index_style, objectiveFuncWrapper, constraintFuncWrapper, 
                           gradientFuncWrapper, constraintJacobianWrapper, hessianWrapper);

        // If you want to verify which options are getting set in the optimizer, you can create a file ipopt.opt
        // with "print_user_options yes", and set print_level to (at least 1).  It will then print the options to the screen.
        AddIpoptNumOption(nlp, "tol", convergenceTolerance);
        AddIpoptIntOption(nlp, "max_iter", maxIterations);
        AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
        AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory"); // needs to be limited-memory unless you have explicit hessians
        AddIpoptIntOption(nlp, "limited_memory_max_history", limitedMemoryHistory);
        AddIpoptIntOption(nlp, "print_level", diagnosticsLevel); // default is 4

        static const char *advancedRealOptions[] = {"obj_scaling_factor", "nlp_scaling_max_gradient",0}; 
        Real value;
        for(int i=0;advancedRealOptions[i];i++) {
            if(getAdvancedRealOption(advancedRealOptions[i],value))
                AddIpoptNumOption(nlp, advancedRealOptions[i], value);
        }

        // Only makes sense to do a warm start if this is not the first call to optimize() (since we need 
        // reasonable starting multiplier values)
        bool use_warm_start=false;
        if(getAdvancedBoolOption("warm_start", use_warm_start) && use_warm_start && !firstOptimization) {
            AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
            AddIpoptStrOption(nlp, "warm_start_entire_iterate", "yes");
            //AddIpoptStrOption(nlp, "warm_start_same_structure", "yes"); // couldn't get this one to work
        } 

        double obj;

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
            char buf[1024];
            sprintf(buf, "Ipopt: %s (status %d)",applicationReturnStatusToString(status).c_str(),status);
            SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf));
        }

        firstOptimization = false;

        return(obj);
    }


} // namespace SimTK
