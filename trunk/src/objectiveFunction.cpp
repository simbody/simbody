#include "objectiveFunction.h"

 namespace SimTK {

// virtual MeasureType GetValue( const ParametersType & parameters ) const = 0;

   virtual double GetValue( double *parameters ) const = 0;

// virtual void GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const = 0;

  virtual void GetGradient( double *parameters, double *gradient ) const = 0;

  virtual void GetValueAndGradient( double *parameters, double *value, double *gradient,) const  
  {
    *value = this->GetValue( parameters );
    this->GetGradient( parameters, gradient );
  };

} // namespace SimTK


