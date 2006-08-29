
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
#include <iostream>
#include "Simmath.h"
#include "Optimizer.h"
#include "OptimizerImplementation.h"

using std::cout;
using std::endl;
extern void lbfgs_( int*n, int*m, double *x, double *f, double *g,
    int *diagco, double *diag, int *iprint, double *eps,
    double *xtol, double *w, int *iflag);

const int NUMBER_OF_CORRECTIONS = 5;

namespace SimTK {

      OptimizerImplementation::OptimizerImplementation( int n ) { // constructor
          smStatus status;
          int m;
      /* internal flags for LBFGS */
         iprint[0] = 1; iprint[1] = 0; // no output generated
         xtol[0] = 1e-16; // from itk/core/vnl/algo/vnl_lbfgs.cxx
         diagco[0] = 0;      // do not supply diagonal of hessian
         user_data = 0;
         objFunc   = 0;

// TODO THROW bad_alloc exception if n<1 || malloc fails
          dimension = n;
          numCorrections = m = NUMBER_OF_CORRECTIONS;
          Algorithm = LBFGS;
          work = new double[dimension*(2*m+1) + 2*m];
          diag = new double[dimension];
          gradient = new SimTK::Vector(dimension);
         
      }

      OptimizerImplementation::OptimizerImplementation() { // constructor
         dimension = 0;
         user_data = 0;
         objFunc   = 0;
      }

      unsigned int OptimizerImplementation::optParamStringToValue( char *parameter )  {

         unsigned int param;

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
           param = UNKNOWN_PARAMETER;
         }

         return( param );

      }

// TODO set pointer to user_data

     smStatus OptimizerImplementation::setOptimizerParameters(unsigned  int parameter, double *values) {
          unsigned int status = SUCCESS;
          int i;

          switch( parameter) {
             case TRACE:
                   Trace = (unsigned int)values[0];
                   break;
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
                   status = UNKNOWN_PARAMETER;
                   break;
          }

  
       return(status); 

      }
      smStatus OptimizerImplementation::getOptimizerParameters(unsigned int parameter, double *values) {
          int status = SUCCESS;
          int i;

            switch( parameter) {
               case TRACE:
                     values[0] = (double )Trace;
                     break;
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
                     status = UNKNOWN_PARAMETER;
                     break;
            }

  
          return(status); 
      }

      smStatus OptimizerImplementation::setObjectiveFunction(SimTK::ObjectiveFunction *objectiveFunction ) {
            objFunc = objectiveFunction;
            return(SUCCESS);
      }

      smStatus OptimizerImplementation::setObjectiveFunction( void (*func)(int, double*,double*,double*, void*)) {

          costFunction = func;
          return(SUCCESS);
      }
     smStatus OptimizerImplementation::optimize( double *results ) {

         smStatus  status = SUCCESS;
         int i;
         int run_optimizer = 1;
         int iflag[1] = {0};
         double f;
         SimTK::Vector &gradient_ref =  *gradient;

         while( run_optimizer ) {

            (*costFunction)( dimension, results, &f, &gradient_ref[0], user_data );

            lbfgs_( &dimension, &numCorrections, &results[0], &f, &gradient_ref[0],
            diagco, diag, iprint, &GradientConvergenceTolerance,  xtol,
            work, iflag );


            if( iflag[0] <= 0 ) {
              run_optimizer = 0;
              status = iflag[0];
            }
         }



         return(status);

      }


     smStatus OptimizerImplementation::optimize(  SimTK::Vector &results ) {

         smStatus  status = SUCCESS;
         int i;
         int run_optimizer = 1;
         int iflag[1] = {0};
         double f;
         SimTK::Vector &gradient_ref =  *gradient;

         while( run_optimizer ) {

            f = objFunc->getValueAndGradient( results, gradient_ref );
            lbfgs_( &dimension, &numCorrections, &results[0], &f, &gradient_ref[0],
            diagco, diag, iprint, &GradientConvergenceTolerance,  xtol,
            work, iflag );


            if( iflag[0] <= 0 ) {
              run_optimizer = 0;
              status = iflag[0];
            }
         }

         return(status);

      }


} // namespace SimTK

/* instantiate */



