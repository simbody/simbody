#ifndef _SimTK_OPTIMIZER_IMPLEMENTATION_H_
#define _SimTK_OPTIMIZER_IMPLEMENTATION_H_

#include <iostream>
#include "Simmath.h"
#include "optimizerInterface.h"
using std::cout;
using std::endl;
extern void lbfgs_( int*n, int*m, double *x, double *f, double *g,
    int *diagco, double *diag, int *iprint, double *eps,
    double *xtol, double *w, int *iflag);

namespace SimTK {

class optimizerImplementation : public smOptimizerInterface {
    public:
      optimizerImplementation( int );
      optimizerImplementation();
      unsigned int optParamStringToValue( char *parameter );
      ~optimizerImplementation(){
          if(dimension > 0 ) {
             delete [] work;
             delete [] g;
             delete [] diag;
             delete [] results;
         }
     }


     smStatus setOptimizerParameters(unsigned  int parameter, double *values) {
          unsigned int status = SUCCESS;
          int i;

          switch( parameter) {
             case INITIAL_VALUES:
                   for(i=0;i<dimension;i++) {
                      results[i] = values[i];
                   }
                   break;
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
/*
     smStatus setOptimizerParameters(unsigned  int param, int *values) {
          cout << "setOptimizerParameters param= " << param << "  values = "<< *values << endl;
          return( smSetOptimizerParametersi ( handle, param, values ));
      }
*/
      smStatus getOptimizerParameters(unsigned int parameter, double *values) {
          int status = SUCCESS;
          int i;

            switch( parameter) {
               case INITIAL_VALUES:
                     for(i=0;i<dimension;i++) {
                        values[i] = results[i];
                     }
                     break;
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
/*
      int getOptimizerParameters(unsigned int param, int *values) {
          cout << "getOptimizerParameters param= " << param << "values = "<< *values << endl;
          return( smGetOptimizerParametersi( handle, param, values ));
      }
*/

      smStatus setObjectiveFunction( void (*func)(double*,double*,double*)) {

          costFunction = func;
          return(SUCCESS);
      }
      smStatus optimize(double *values ) {


         smStatus  status = SUCCESS;
         int i;
         int run_optimizer = 1;
         int iprint[2] = {1, 0}; // no output generated
         int iflag[1] = {0};
         double xtol[1] = {1e-16}; // from itk/core/vnl/algo/vnl_lbfgs.cxx
         int diagco[1] = {0};      // do not supply diagonal of hessian


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
         if( status == SUCCESS ) {
              for( i=0; i<dimension; i++ ) {
              *(values++) = results[i];
              }
         }

         return(status);

      }

  private:

     int         dimension;  // dimension of the problem
     int         numCorrections;
     int         Trace;
     int         Algorithm;
     int         MaxNumFuncEvals;
     double               *work;
     double               *results;
     double               f;
     double               *g;
     double               *diag;
     double               GradientConvergenceTolerance;
     double               LineSearchAccuracy;
     double               DefaultStepLength;
     void                 (*costFunction)(double*, double*, double*);

};
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_IMPLEMENTATION_H_
