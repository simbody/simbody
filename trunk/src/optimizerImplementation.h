#ifndef _SimTK_OPTIMIZER_IMPLEMENTATION_H_
#define _SimTK_OPTIMIZER_IMPLEMENTATION_H_

#include <iostream>
#include "Simmath.h"
#include "optimizerInterface.h"
#include "objectiveFunction.h"

namespace SimTK {

class optimizerImplementation : public smOptimizerInterface {
    public:
     optimizerImplementation( int );
     optimizerImplementation();
     unsigned int optParamStringToValue( char *parameter );
     smStatus setOptimizerParameters(unsigned  int parameter, double *values); 
     smStatus getOptimizerParameters(unsigned int parameter, double *values); 
     smStatus setObjectiveFunction( void (*func)(int,double*,double*,double*,void*)) ;
     smStatus setObjectiveFunction( SimTK::objectiveFunction *objFunc);

     smStatus optimize( double *results ); 

     smStatus optimize(  SimTK::Vector_<Real> &results ); 


      ~optimizerImplementation(){
          if(dimension > 0 ) {
             delete [] work;
             delete [] diag;
         }
     }

     int         dimension;  // dimension of the problem
     int         numCorrections;
     int         Trace;
     int         Algorithm;
     int         MaxNumFuncEvals;
     double     *work;
     double     *diag;
     double      GradientConvergenceTolerance;
     double      LineSearchAccuracy;
     double      DefaultStepLength;
     void        (*costFunction)(int, double*, double*, double*, void*);
     int         iprint[2]; 
     double      xtol[1]; 
     int         diagco[1];
     void       *user_data;
     SimTK::Vector_<Real> *gradient;
     SimTK::objectiveFunction *objFunc;


};
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_IMPLEMENTATION_H_
