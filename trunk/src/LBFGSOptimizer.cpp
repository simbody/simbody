

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
    double *xtol, double *w, int *iflag, int *info);

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

// TODO option for xtol  ??

         if( sys.getNumParameters() < 1 ) {
             char *where = "Optimizer Initialization";
             char *szName= "dimension";
             SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  sys.getNumParameters(), INT_MAX, where); 
         }
         n = sys.getNumParameters();
         m = NUMBER_OF_CORRECTIONS;
         work = new double[n*(2*m+1) + 2*m];
         diag = new double[n];
         gradient = new double[n];
     } 


     double LBFGSOptimizer::optimize(  Vector &results ) {

         int i,info;
         int iflag[1] = {0};
         double f;
         const OptimizerSystem& sys = getOptimizerSystem();
         int n = sys.getNumParameters();
         int m = NUMBER_OF_CORRECTIONS;
         char buf[256];

         iprint[1] = diagnosticsLevel;


         do {   

            objectiveFuncWrapper( n, &results[0], true, &f, (void*)this );
            gradientFuncWrapper( n,  &results[0], false, gradient, (void*)this );

            lbfgs_( &n, &m, &results[0], &f, gradient,
            diagco, diag, iprint, &convergenceTolerance,  xtol,
            work, iflag, &info );

         } while( iflag[0] > 0 );

         if( iflag[0] == 0 ) {
            objectiveFuncWrapper( n, &results[0], true, &f, (void*)this );
            return f;
         } else {
            if( iflag[0] == -1 ) {
               if( info == 0) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: IMPROPER INPUT PARAMETERS");
               } else if( info == 2) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL");
               } else if( info == 3) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: MORE THAN 20 FUNCTION EVALUATIONS WERE REQUIRED AT THE PRESENT ITERATION");
               } else if( info == 4) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: THE STEP IS TOO SMALL");
               } else if( info == 5) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: THE STEP IS TOO LARGE");
               } else if( info == 6){
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: ROUNDING ERRORS PREVENT FURTHER PROGRESS.\n THERE MAY NOT BE A STEP WHICH SATISFIES THE SUFFICIENT DECREASE AND CURVATURE\n CONDITIONS. TOLERANCES MAY BE TOO SMALL.");
               } else {
                  sprintf(buf, "LBFGS ERROR: info = %d \n",info );
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf) );
               }
            } else if ( iflag[0] == -2 ) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: The i-th diagonal element of the diagonal inverse Hessian approximation, \n given in DIAG, is not positive.");
            } else if ( iflag[0] == -3 ) {
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, 
                  "LBFGS ERROR: Improper input parameters for LBFGS (N or M are ot positive). \n" );
            } else {
                  sprintf(buf, "LBFGS ERROR: iflag = %d \n",iflag[0] );
                  SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf) );
            }
      }

   }
} // namespace SimTK
