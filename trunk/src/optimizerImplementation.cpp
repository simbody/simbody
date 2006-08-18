#include <iostream>
#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

using std::cout;
using std::endl;
const int NUMBER_OF_CORRECTIONS = 5;

namespace SimTK {

      optimizerImplementation::optimizerImplementation( int n ) { // constructor
          smStatus status;
          int m;
// TODO THROW bad_alloc exception if n<1 || malloc fails
          dimension = n;
          numCorrections = m = NUMBER_OF_CORRECTIONS;
          Algorithm = LBFGS;
          work = new double[dimension*(2*m+1) + 2*m];
          g =    new double[dimension];
          diag = new double[dimension];
          results = new double[dimension];
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
         } else if( 0 == strncmp( "INITIAL_VALUES", parameter, 1)) {
           param = INITIAL_VALUES;
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


} // namespace SimTK

