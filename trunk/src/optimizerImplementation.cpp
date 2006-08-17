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


} // namespace SimTK

