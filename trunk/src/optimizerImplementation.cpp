#include <iostream>
#include "Simmath.h"
#include "optimizer.h"
#include "optimizerImplementation.h"

using std::cout;
using std::endl;

namespace SimTK {

      optimizerImplementation::optimizerImplementation( int n ) { // constructor
          dim = n;
          cout << "optimizerImplementation constructor n = " <<  dim << endl;
      }

} // namespace SimTK

