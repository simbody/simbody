//-----------------------------------------------------------------------------
// File:     UserSpringDampingForces.h
// Class:    UserSpringDampingForces
// Parent:   GeneralForceElements
// Children: None
// Purpose:  Applies Spring-Damper forces to the model.
// Author:   Ajay Seth - May 31, 2007
//-----------------------------------------------------------------------------

#include "SimTKsimbody.h"
using namespace SimTK;

//--------------------------------------------------------------------------
// User-defined classes for adding forces/torques are constructed as follows:
// 1. Create a constructor with whatever arguments make sense for the force or torque and copy the arguments into class data.
// 2. Create a clone method (all clone methods are identical except for the class name appearing after "new")
// 3. Create a calc  method (the arguments and return type for all calc methods are identical).
//    The code in the calc method is specific to the calculation of force or torque.
// Note: The set of all forces is replaced by an equivalent set, consisting of a torque
//       that is equal to the moment of the forces about the body's origin together
//       with the resultant of the forces applied at the body's origin.
//--------------------------------------------------------------------------
class UserSpringDampingForces : public GeneralForceElements::UserForce 
{
public:
   // Constructor is explicit
   explicit UserSpringDampingForces(Real stiffness, Real damping)
   { 
      k = stiffness;
      d = damping;
   } 

   // The clone method is used internally by Simbody (required by virtual parent class)
   UserForce*  clone() const  { return new UserSpringDampingForces(*this); }

   // The calc method is where forces or torques are calculated (required by virtual parent class)
   void  calc( const MatterSubsystem& matter,              // Input information (matter)
               const State&           state,               // Input information (current state)
               Vector_<SpatialVec>&   bodyForces,          // Forces and torques on bodies
               Vector_<Vec3>&         particleForces,      // Forces on particles (currently unused)
               Vector&                mobilityForces,      // Generalized forces 
               Real&                  pe ) const           // For forces with a potential energy
    {
        // Query the matter subsystem for the joint angles
		const Vector&  q  = state.getQ();
		const Vector&  u  = state.getU();

		mobilityForces += (-k*q-d*u);
		pe += 0.5*k*(~q*q);
    }

private:
    Real k;
    Real d;
};