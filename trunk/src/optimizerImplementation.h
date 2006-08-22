#ifndef _SimTK_OPTIMIZER_IMPLEMENTATION_H_
#define _SimTK_OPTIMIZER_IMPLEMENTATION_H_

#include <iostream>
#include "Simmath.h"
#include "optimizerInterface.h"

namespace SimTK {

class optimizerImplementation : public smOptimizerInterface {
    public:
     optimizerImplementation( int );
     optimizerImplementation();
     unsigned int optParamStringToValue( char *parameter );
     smStatus setOptimizerParameters(unsigned  int parameter, double *values); 
     smStatus getOptimizerParameters(unsigned int parameter, double *values); 
     smStatus setObjectiveFunction( void (*func)(double*,double*,double*)) ;

     template < typename T >
     smStatus optimize( T results ); 

     template < int N, typename T, int S >
     smStatus optimize(  SimTK::Vec< N, T, S> results ); 


      ~optimizerImplementation(){
          if(dimension > 0 ) {
             delete [] work;
             delete [] g;
             delete [] diag;
         }
     }

     int         dimension;  // dimension of the problem
     int         numCorrections;
     int         Trace;
     int         Algorithm;
     int         MaxNumFuncEvals;
     double               *work;
     double               f;
     double               *g;
     double               *diag;
     double               GradientConvergenceTolerance;
     double               LineSearchAccuracy;
     double               DefaultStepLength;
     void                 (*costFunction)(double*, double*, double*);
     int iprint[2]; 
     double xtol[1]; 
     int diagco[1];


};
} // namespace SimTK
#endif //_SimTK_OPTIMIZER_IMPLEMENTATION_H_
