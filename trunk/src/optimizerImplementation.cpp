#include <iostream>
#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

using std::cout;
using std::endl;
extern void lbfgs_( int*n, int*m, double *x, double *f, double *g,
    int *diagco, double *diag, int *iprint, double *eps,
    double *xtol, double *w, int *iflag);

const int NUMBER_OF_CORRECTIONS = 5;

namespace SimTK {

      optimizerImplementation::optimizerImplementation( int n ) { // constructor
          smStatus status;
          int m;
         iprint[0] = 1; iprint[1] = 0; // no output generated
         xtol[0] = 1e-16; // from itk/core/vnl/algo/vnl_lbfgs.cxx
         diagco[0] = 0;      // do not supply diagonal of hessian


// TODO THROW bad_alloc exception if n<1 || malloc fails
          dimension = n;
          numCorrections = m = NUMBER_OF_CORRECTIONS;
          Algorithm = LBFGS;
          work = new double[dimension*(2*m+1) + 2*m];
          g =    new double[dimension];
          diag = new double[dimension];
      }

      optimizerImplementation::optimizerImplementation( ) { // constructor
         dimension = 0;
      }

      unsigned int optimizerImplementation::optParamStringToValue( char *parameter )  {

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



     smStatus optimizerImplementation::setOptimizerParameters(unsigned  int parameter, double *values) {
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
      smStatus optimizerImplementation::getOptimizerParameters(unsigned int parameter, double *values) {
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

      smStatus optimizerImplementation::setObjectiveFunction( void (*func)(double*,double*,double*)) {

          costFunction = func;
          return(SUCCESS);
      }
     template < typename T >
     smStatus optimizerImplementation::optimize( T results ) {

         smStatus  status = SUCCESS;
         int i;
         int run_optimizer = 1;
         int iflag[1] = {0};


         while( run_optimizer ) {

            (*costFunction)( results, &f, g );

            lbfgs_( &dimension, &numCorrections, results, &f, g,
            diagco, diag, iprint, &GradientConvergenceTolerance,  xtol,
            work, iflag );


            if( iflag[0] <= 0 ) {
              run_optimizer = 0;
              status = iflag[0];
            }
         }

         return(status);

      }


    template < int N, typename T, int S >
     smStatus optimizerImplementation::optimize(  SimTK::Vec< N, T, S> results ) {

         smStatus  status = SUCCESS;
         int i;
         int run_optimizer = 1;
         int iflag[1] = {0};


         while( run_optimizer ) {

            (*costFunction)( &results[0], &f, g );

            lbfgs_( &dimension, &numCorrections, &results[0], &f, g,
            diagco, diag, iprint, &GradientConvergenceTolerance,  xtol,
            work, iflag );


            if( iflag[0] <= 0 ) {
              run_optimizer = 0;
              status = iflag[0];
            }
         }

         return(status);

      }

template smStatus optimizerImplementation::optimize(double *);

template smStatus optimizerImplementation::optimize<2, double, 1>(Vec<2, double, 1>);
template smStatus optimizerImplementation::optimize(Vec<3, double, 1>);

} // namespace SimTK

/* instantiate */



