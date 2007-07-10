

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
#include "InteriorPointOptimizer.h"

using std::cout;
using std::endl;


namespace SimTK {


// Assume by the time this constructor is called, the number of parameters and constraints has been finalized
InteriorPointOptimizer::InteriorPointOptimizer( OptimizerSystem& sys )
        : OptimizerRep( sys ) {

        int n = sys.getNumParameters();
        int m = sys.getNumConstraints();

        if( n < 1 ) {
            char where[] = " InteriorPointOptimizer Initialization";
            char szName[] = "dimension";
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
        char tol[] = "tol";
        AddIpoptNumOption(nlp, tol, convergenceTolerance);
        char max_iter[] = "max_iter";
        AddIpoptIntOption(nlp, max_iter, maxIterations);
        char mu_strategy[] = "mu_strategy";
        char adaptive[] = "adaptive";
        AddIpoptStrOption(nlp, mu_strategy, adaptive);
        char hessian_approximation[] = "hessian_approximation";
        char limited_memory[] = "limited-memory";
        AddIpoptStrOption(nlp, hessian_approximation, limited_memory); // needs to be limited-memory unless you have explicit hessians
        char limited_memory_max_history[] = "limited_memory_max_history";
        AddIpoptIntOption(nlp, limited_memory_max_history, limitedMemoryHistory);
        char print_level[] = "print_level";
        AddIpoptIntOption(nlp, print_level, diagnosticsLevel); // default is 4

        // should be const char * but AddIpoptNumOption expects a char *
        static char *advancedRealOptions[3];
        char obj_scaling_factor[] = "obj_scaling_factor";
        advancedRealOptions[0] = obj_scaling_factor;
        char nlp_scaling_max_gradient[] = "nlp_scaling_max_gradient";
        advancedRealOptions[1] = nlp_scaling_max_gradient;
        advancedRealOptions[2] = 0;
       
        Real value;
        for(int i=0;advancedRealOptions[i];i++)
            if(getAdvancedRealOption(advancedRealOptions[i],value))
                AddIpoptNumOption(nlp, advancedRealOptions[i], value);

        // Only makes sense to do a warm start if this is not the first call to optimize() (since we need 
        // reasonable starting multiplier values)
        bool use_warm_start=false;
        char warm_start[] = "warm_start";
        if(getAdvancedBoolOption(warm_start, use_warm_start) && use_warm_start && !firstOptimization) {
            char warm_start_init_point[] = "warm_start_init_point";
            char warm_start_entire_iterate[] =  "warm_start_entire_iterate"; 
            char yes_string[] = "yes";
            AddIpoptStrOption(nlp, warm_start_init_point, yes_string);
            AddIpoptStrOption(nlp, warm_start_entire_iterate, yes_string);
            //AddIpoptStrOption(nlp, "warm_start_same_structure", "yes"); // couldn't get this one to work
        } 

        double obj;

        int status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, (void *)this );

        if (status != Solve_Succeeded) {
            char buf[1024];
            sprintf(buf, "Ipopt failed with status = %d",status );
            SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf));
        }

        FreeIpoptProblem(nlp); 

        firstOptimization = false;

        // Only delete these if they aren't pointing to existing parameter limits
        if( !getOptimizerSystem().getHasLimits() ) {
           delete [] x_U;
           delete [] x_L;
        }

        return(obj);
    }


} // namespace SimTK
