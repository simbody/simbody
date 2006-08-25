#include "objectiveFunction.h"

 namespace SimTK {


   double getValue( SimTK::Vector_<Real>&) const = 0;


  void getGradient(  SimTK::Vector_<Real>&, SimTK::Vector_<Real>&) const = 0;

  double getValueAndGradient( SimTK::Vector_<Real>& parameters,  SimTK::Vector_<Real>& gradient,) const  
  {
    this->GetGradient( parameters, gradient );
    return( this->GetValue( parameters ));
  };

} // namespace SimTK


