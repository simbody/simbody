

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


InteriorPointOptimizer::InteriorPointOptimizer( OptimizerSystem& sys )
        : OptimizerRep( sys ) {

        printf(" InteriorPointOptimizer constructor \n");         
          int i,n;
          char buf[1024];

          /* C-style; start counting of rows and column indices at 0 */
          Index index_style = 0; 


         if( sys.dimension < 1 ) {
             char *where = " InteriorPointOptimizer Initialization";
             char *szName= "dimension";
             SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  sys.dimension, INT_MAX, where); 
         } else {
            n = sys.dimension;
            /* set the bounds on the equality constraint functions */
            m = sys.numConstraints;
            g_U = (double *)malloc(sizeof(double)*m); // TODO free these
            g_L = (double *)malloc(sizeof(double)*m);
            double *x_U = sys.upper_bounds;
            double *x_L = sys.lower_bounds;
            for(i=0;i<sys.numEqualConstraints;i++){
                g_U[i] = g_L[i] = 0.0;
            }

            Index nele_hess = 0;
            Index nele_jac = n*m; /* always assume dense
            
            /* set the bounds on the inequality constraint functions */
            for(i=sys.numEqualConstraints;i<m;i++){
                g_U[i] = 2e19;
                g_L[i] = 0.0;
            }

            mult_x_L = (Number*)malloc(sizeof(Number)*n);
            mult_x_U = (Number*)malloc(sizeof(Number)*n);
        
            nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, 
                  nele_hess, index_style, objectiveFuncWrapper, constraintFuncWrapper, 
                  gradientFuncWrapper, constraintJacobianWrapper, hessianWrapper);

            AddIpoptNumOption(nlp, "tol", 1e-3);
            AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
            AddIpoptStrOption(nlp, "output_file", "ipopt.out");
            AddIpoptStrOption(nlp, "linear_solver", "lapack");
            AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
            AddIpoptIntOption(nlp, "print_level", 4); // default is 4

          }


     } 


     double InteriorPointOptimizer::optimize(  Vector &results ) {

         int i;
         double obj;
         double *x = &results[0];
         int status;

        printf("call InteriorPoint optimize \n");

         status = IpoptSolve(nlp, x, NULL, &obj, NULL, mult_x_L, mult_x_U, (void *)this );

         if (status == Solve_Succeeded) {
             printf("Ipopt CONVERGED \n");
         } else {
             printf("Ipopt Solve failed  \n");
         }

         return(obj);
      }

      unsigned int InteriorPointOptimizer::optParamStringToValue( char *parameter )  {

         unsigned int param;
         char buf[1024];

         if( 0 == strncmp( "FUNCION_EVALUATIONS", parameter, 1) ) {
           param = MAX_FUNCTION_EVALUATIONS;
         } else if( 0 == strncmp( "STEP_LENGTH", parameter, 1)) {
           param = DEFAULT_STEP_LENGTH;
         } else if( 0 == strncmp( "TOLERANCE", parameter, 1)) {
           param = TRACE;
         } else if( 0 == strncmp( "GRADIENT", parameter, 1)) {
           param = GRADIENT_CONVERGENCE_TOLERANCE;
         } else if( 0 == strncmp( "ACCURACY", parameter, 1)) {
           param = LINE_SEARCH_ACCURACY;
         } else {
             sprintf(buf," Parameter=%s",parameter);
             SimTK_THROW1(SimTK::Exception::UnrecognizedParameter, SimTK::String(buf) ); 
         }

         return( param );

      }

     void InteriorPointOptimizer::setOptimizerParameters(unsigned int parameter, double *values ) { 
          int i;
          char buf[1024];

          printf("call InteriorPoint setOptimizerParameters \n");
          switch( parameter) {
             case MAX_FUNCTION_EVALUATIONS:
                   MaxNumFuncEvals = (unsigned int)values[0];
                   break;
             case DEFAULT_STEP_LENGTH:
                   DefaultStepLength = (unsigned int)values[0];
                   break;
             case LINE_SEARCH_ACCURACY:
                   LineSearchAccuracy = values[0];
                   break;
             case  GRADIENT_CONVERGENCE_TOLERANCE:
                   GradientConvergenceTolerance = values[0];
                   break;
             default:
/*  TODO fix this 
                   sprintf(buf," Parameter=%d",parameter);
                   SimTK_THROW1(SimTK::Exception::UnrecognizedParameter, SimTK::String(buf) ); 
*/
                   break;
          }

        return; 

      }
     void InteriorPointOptimizer::getOptimizerParameters(unsigned int parameter, double *values ) {
          int i;
          char buf[1024];

        printf("call InteriorPoint getOptimizerParameters \n");

            switch( parameter) {
               case MAX_FUNCTION_EVALUATIONS:
                     values[0] = (double )MaxNumFuncEvals;
                     break;
               case DEFAULT_STEP_LENGTH:
                     values[0] = DefaultStepLength;
                     break;
               case LINE_SEARCH_ACCURACY:
                     values[0] = LineSearchAccuracy;
                     break;
               case  GRADIENT_CONVERGENCE_TOLERANCE:
                      values[0] = GradientConvergenceTolerance;
                      break;
               default:
                      sprintf(buf," Parameter=%d",parameter);
                      SimTK_THROW1(SimTK::Exception::UnrecognizedParameter, SimTK::String(buf) ); 
                      break;
            }
  
          return; 
   }
} // namespace SimTK
