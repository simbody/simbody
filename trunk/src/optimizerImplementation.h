#ifndef _SimTK_OPTIMIZER_IMPLEMENTATION_H_
#define _SimTK_OPTIMIZER_IMPLEMENTATION_H_

#include <iostream>
#include "Simmath.h"
#include "optimizerInterface.h"
using std::cout;
using std::endl;

namespace SimTK {

class optimizerImplementation : public smOptimizerInterface {
    public:
      optimizerImplementation( int );
      optimizerImplementation();
      
      int dim;  // dimension of the problem

     int setOptimizerParameters(int param, double *values) {
          cout << "setOptimizerParameters param= " << param << "  values = "<< *values << endl;
          return(SUCCESS);
      }

      int getOptimizerParameters(int param, double *values) {
          cout << "getOptimizerParameters param= " << param << "values = "<< *values << endl;
          return(SUCCESS);
      }

      int setObjectiveFunction(int *func ) {
          cout << "setCostFunction = " <<  endl;
          return(SUCCESS);
      }
      int setInitialPosition(double *pos ) {
          cout << "setInitialPosition = " << *pos <<  endl;
          return(SUCCESS);
      }
      int optimize(double *results ) {
          results[0] = 5.67;
          cout << "optimize = " << *results <<  endl;
          return(SUCCESS);
      }


};
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_IMPLEMENTATION_H_
