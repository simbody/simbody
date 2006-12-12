

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
#include "LBFGSOptimizer.h"

using std::cout;
using std::endl;
extern void lbfgs_( int*n, int*m, double *x, double *f, double *g,
    int *diagco, double *diag, int *iprint, double *eps,
    double *xtol, double *w, int *iflag);

namespace SimTK {
const int NUMBER_OF_CORRECTIONS = 5;

     LBFGSOptimizer::LBFGSOptimizer( OptimizerSystem& sys )
        : OptimizerRep( sys ) 
{
          int n,m;
          char buf[1024];


      /* internal flags for LBFGS */
         iprint[0] = 1; iprint[1] = 0; // no output generated
         xtol[0] = 1e-16; // from itk/core/vnl/algo/vnl_lbfgs.cxx
         diagco[0] = 0;      // do not supply diagonal of hessian

         if( sys.getNumParameters() < 1 ) {
             char *where = "Optimizer Initialization";
             char *szName= "dimension";
             SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  sys.getNumParameters(), INT_MAX, where); 
         }
         n = sys.getNumParameters();
         numCorrections = m = NUMBER_OF_CORRECTIONS;
         work = new double[n*(2*m+1) + 2*m];
         diag = new double[n];
         gradient = new double[n];
     } 


     double LBFGSOptimizer::optimize(  Vector &results ) {

         int i;
         int run_optimizer = 1;
         int iflag[1] = {0};
         double f;
         const OptimizerSystem& sys = getOptimizerSystem();
         int n = sys.getNumParameters();

         while( run_optimizer ) {  // TODO callbacks use ptr to functions  

            objectiveFuncWrapper( n, &results[0], true, &f, (void*)this );
            gradientFuncWrapper( n,  &results[0], false, gradient, (void*)this );

            lbfgs_( &n, &numCorrections, &results[0], &f, gradient,
            diagco, diag, iprint, &GradientConvergenceTolerance,  xtol,
            work, iflag );

            /* TODO check iflag[0] for status errors */
            if( iflag[0] <= 0 ) {
              run_optimizer = 0;
            }
         }
         objectiveFuncWrapper( n, &results[0], true, &f, (void*)this );
         return f;
      }

      unsigned int LBFGSOptimizer::optParamStringToValue( char *parameter )  {

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

     void LBFGSOptimizer::setOptimizerParameters(unsigned int parameter, double *values ) { 
          int i;
          char buf[1024];

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
                   sprintf(buf," Parameter=%d",parameter);
                   SimTK_THROW1(SimTK::Exception::UnrecognizedParameter, SimTK::String(buf) ); 
                   break;
          }

        return; 

      }
     void LBFGSOptimizer::getOptimizerParameters(unsigned int parameter, double *values ) {
          int i;
          char buf[1024];


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
                      sprintf(buf," Parameter=%d",parameter);
                      SimTK_THROW1(SimTK::Exception::UnrecognizedParameter, SimTK::String(buf) ); 
                      break;
            }
  
          return; 
   }
} // namespace SimTK
