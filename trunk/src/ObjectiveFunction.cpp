#include "ObjectiveFunction.h"

 namespace SimTK {


  double getValue( SimTK::Vector&) const = 0;


  void getGradient(  SimTK::Vector&, SimTK::Vector&) const = 0;

  double getValueAndGradient( SimTK::Vector& parameters,  SimTK::Vector& gradient,) const  
  {
    this->GetGradient( parameters, gradient );
    return( this->GetValue( parameters ));
  };

} // namespace SimTK


