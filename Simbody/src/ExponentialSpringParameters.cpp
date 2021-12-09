/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021 Authors.
 Authors: Frank C. Anderson
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ----------------------------------------------------------------------------*/


#include "SimTKcommon.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/ExponentialSpringForce.h"

using namespace SimTK;
using std::cout;
using std::endl;

//=============================================================================
// Class - ExponentialSpringParameters
//=============================================================================
//_____________________________________________________________________________
// Default Constructor
ExponentialSpringParameters::
ExponentialSpringParameters() :
    d0(0.0065905), d1(0.5336), d2(1150.0), kvNorm(0.5),
    kpFric(20000.0), kvFric(0.0),
    tau(0.01), vSettle(0.001), aSettle(0.1) {

    setElasticityAndViscosityForCriticalDamping(kpFric);
}
//_____________________________________________________________________________
// Copy Constructor
ExponentialSpringParameters::
ExponentialSpringParameters(const ExponentialSpringParameters& source) {
    operator=(source);
}
//_____________________________________________________________________________
// Copy Assignment Operator
ExponentialSpringParameters&
ExponentialSpringParameters::
operator=(const ExponentialSpringParameters& source) {
    if(&source != this) {
        d0 = source.d0;
        d1 = source.d1;
        d2 = source.d2;
        kvNorm = source.kvNorm;
        kpFric = source.kpFric;
        kvFric = source.kvFric;
        tau = source.tau;
        vSettle = source.vSettle;
        aSettle = source.aSettle;
    }
    return *this;
}
//_____________________________________________________________________________
// Equality Operator
bool
ExponentialSpringParameters::
operator==(const ExponentialSpringParameters& other) const {
    if(&other == this) return true;
    return ((d0 == other.d0) && (d1 == other.d1) && (d2 == other.d2) &&
        (kvNorm == other.kvNorm) && (kpFric == other.kpFric) &&
        (tau == other.tau) &&
        (vSettle == other.vSettle) && (aSettle == other.aSettle));
}
//_____________________________________________________________________________
// Set the parameters that control the shape of the exponential function.
void
ExponentialSpringParameters::
setShapeParameters(Real d0, Real d1, Real d2) {
    // d0
    this->d0 = d0;

    // d1
    if(d1 <= 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - d1 should be positive!"
            << endl;
    } else this->d1 = d1;

    // d2
    if(d2 <= 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - d2 should be positive!"
            << endl;
    } else this->d2 = d2;
}
//_____________________________________________________________________________
// Get the parameters that control the shape of the exponential function.
void
ExponentialSpringParameters::
getShapeParameters(Real &d0, Real& d1, Real& d2) const {
    d0 = this->d0;
    d1 = this->d1;
    d2 = this->d2;
}
//_____________________________________________________________________________
// Set the viscosity of the exponential spring.  This quanitiy only affects
// the damping in the direction normal to the floor.
void
ExponentialSpringParameters::
setNormalViscosity(Real kvNorm) {
    if(kvNorm < 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - kvNorm should be zero "
            << "or positive!" << endl;
    } else this->kvNorm = kvNorm;
}
//_____________________________________________________________________________
// Get the viscosity of the exponential spring.  This is the viscosity
// that applies to velocity in the normal direction.
Real
ExponentialSpringParameters::
getNormalViscosity() const {
    return kvNorm;
}
//_____________________________________________________________________________
// Set the elasticity and compute the visosity to produce critically
// damped motion for a specified mass.
void
ExponentialSpringParameters::
setElasticityAndViscosityForCriticalDamping(Real kp, Real mass) {
    // Set the elasticity
    setElasticity(kp);

    // Compute the viscosity
    if(mass<=0.0) {
        // An exception should be throw, but for now...
        cout<<"ExponentialSpringParameters: ERR - mass should be positive!"
            <<endl;
    } else this->kvFric = 2.0 * std::sqrt(this->kpFric * mass);
}
//_____________________________________________________________________________
// Set the elasticity of the friction spring.
void
ExponentialSpringParameters::
setElasticity(Real kp) {
    if(kp <= 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - kpFric should be positive!"
            << endl;
    } else this->kpFric = kp;
}
//_____________________________________________________________________________
// Get the elasticity of the friction spring.
Real
ExponentialSpringParameters::
getElasticity() const {
    return kpFric;
}
//_____________________________________________________________________________
// Set the viscosity of the friction spring.
void
ExponentialSpringParameters::
setViscosity(Real kv) {
    if(kv < 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - kvFric should be zero or "
            << "positive!" << endl;
    } else this->kvFric = kv;
}
//_____________________________________________________________________________
// Get the viscosity of the friction spring.
Real
ExponentialSpringParameters::
getViscosity() const {
    return kvFric;
}
//_____________________________________________________________________________
// Set the time constant for transitioning between kinetic and static
// frictional coefficients.
void
ExponentialSpringParameters::
setSlidingTimeConstant(Real tau) {
    if(tau <= 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - tau should be positive!"
            << endl;
    } else this->tau = tau;
}
//_____________________________________________________________________________
// Get the elasticity of the friction spring.
Real
ExponentialSpringParameters::
getSlidingTimeConstant() const {
    return tau;
}
//_____________________________________________________________________________
// Set the velocity for settling into using the static coefficient of friction.
void
ExponentialSpringParameters::
setSettleVelocity(Real vSettle) {
    if(vSettle<=0.0) {
        // An exception should be throw, but for now...
        cout<<"ExponentialSpringParameters: ERR - vSettle should be positive!"
            <<endl;
    } else this->vSettle = vSettle;
}
//_____________________________________________________________________________
// Get the settle velocity.
Real
ExponentialSpringParameters::
getSettleVelocity() const {
    return vSettle;
}
//_____________________________________________________________________________
// Set the acceleration for settling into using the static coefficient of
// friction.
void
ExponentialSpringParameters::
setSettleAcceleration(Real aSettle) {
    if(aSettle <= 0.0) {
        // An exception should be throw, but for now...
        cout << "ExponentialSpringParameters: ERR - aSettle should be positive!"
            << endl;
    } else this->aSettle = aSettle;
}
//_____________________________________________________________________________
// Get the settle acceleration.
Real
ExponentialSpringParameters::
getSettleAcceleration() const {
    return aSettle;
}

