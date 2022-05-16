/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021-22 Authors.
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
ExponentialSpringParameters::
ExponentialSpringParameters() :
    d0(0.0065905), d1(0.5336), d2(1150.0), cz(0.5), maxFz(100000.0),
    kxy(20000.0), cxy(0.0),
    tau(0.01), vSettle(0.01), aSettle(1.0),
    initMus(0.7), initMuk(0.5) {

    setElasticityAndViscosityForCriticalDamping(kxy);
}
//_____________________________________________________________________________
bool
ExponentialSpringParameters::
operator==(const ExponentialSpringParameters& other) const {
    if(&other == this) return true;
    return ((d0 == other.d0) && (d1 == other.d1) && (d2 == other.d2) &&
        (cz == other.cz) && (maxFz == other.maxFz) &&
        (kxy == other.kxy) && (cxy == other.cxy) &&
        (tau == other.tau) &&
        (vSettle == other.vSettle) && (aSettle == other.aSettle) &&
        (initMus == other.initMus) && (initMuk == other.initMuk));
}
//_____________________________________________________________________________
bool
ExponentialSpringParameters::
operator!=(const ExponentialSpringParameters& other) const {
    return !(*this == other);
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setShapeParameters(Real d0, Real d1, Real d2) {
    // d0
    this->d0 = d0;

    // d1
    SimTK_APIARGCHECK1_ALWAYS(d1 > 0.0,
        "ExponentialSpringParameters", "setShapeParameters",
        "expected d1 > 0.0, but d1 = %lf", d1);
    this->d1 = d1;

    // d2
    SimTK_APIARGCHECK1_ALWAYS(d2 > 0.0,
        "ExponentialSpringParameters", "setShapeParameters",
        "expected d2 > 0.0, but d2 = %lf", d2);
    this->d2 = d2;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
getShapeParameters(Real& d0, Real& d1, Real& d2) const {
    d0 = this->d0;
    d1 = this->d1;
    d2 = this->d2;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setNormalViscosity(Real cz) {
    SimTK_APIARGCHECK1_ALWAYS(cz >= 0.0,
        "ExponentialSpringParameters", "setNormalViscosity",
        "expected cz >= 0.0, but cz = %lf", cz);
    this->cz = cz;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getNormalViscosity() const {
    return cz;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setMaxNormalForce(Real maxFz) {
    SimTK_APIARGCHECK1_ALWAYS(maxFz > 0.0,
        "ExponentialSpringParameters", "setMaxNormalForce",
        "expected maxFz > 0.0, but cz = %lf", maxFz);
    this->maxFz = maxFz;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getMaxNormalForce() const {
    return maxFz;
}

//_____________________________________________________________________________
void
ExponentialSpringParameters::
setElasticityAndViscosityForCriticalDamping(Real kxy, Real mass) {
    // Set the elasticity
    setFrictionElasticity(kxy);

    // Compute the viscosity
    SimTK_APIARGCHECK1_ALWAYS(mass > 0.0,
        "ExponentialSpringParameters",
        "setElasticityAndViscosityForCriticalDamping",
        "expected mass > 0.0, but mass = %lf", mass);
    this->cxy = 2.0 * std::sqrt(this->kxy * mass);
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setFrictionElasticity(Real kxy) {
    SimTK_APIARGCHECK1_ALWAYS(kxy > 0.0,
        "ExponentialSpringParameters",
        "setFrictionElasticity",
        "expected kxy > 0.0, but kxy = %lf", kxy);
    this->kxy = kxy;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getFrictionElasticity() const {
    return kxy;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setFrictionViscosity(Real cxy) {
    SimTK_APIARGCHECK1_ALWAYS(cxy >= 0.0,
        "ExponentialSpringParameters",
        "setFrictionViscosity",
        "expected cxy >= 0.0, but cxy = %lf", cxy);
    this->cxy = cxy;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getFrictionViscosity() const {
    return cxy;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setSlidingTimeConstant(Real tau) {
    SimTK_APIARGCHECK1_ALWAYS(tau > 0.0,
        "ExponentialSpringParameters",
        "setSlidingTimeConstant",
        "expected tau > 0.0, but tau = %lf", tau);
    this->tau = tau;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getSlidingTimeConstant() const {
    return tau;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setSettleVelocity(Real vSettle) {
    SimTK_APIARGCHECK1_ALWAYS(vSettle > 0.0,
        "ExponentialSpringParameters",
        "setSettleVelocity",
        "expected vSettle > 0.0, but vSettle = %lf", vSettle);
    this->vSettle = vSettle;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getSettleVelocity() const {
    return vSettle;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setSettleAcceleration(Real aSettle) {
    SimTK_APIARGCHECK1_ALWAYS(aSettle > 0.0,
        "ExponentialSpringParameters",
        "setSettleAcceleration",
        "expected aSettle > 0.0, but aSettle = %lf", aSettle);
    this->aSettle = aSettle;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getSettleAcceleration() const {
    return aSettle;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setInitialMuStatic(Real mus) {
    initMus = mus;
    // Keep initMus greater than or equal to 0.0.
    if(initMus < 0.0) initMus = 0.0;
    // Make sure initMuk is less than or equal to initMus
    if(initMuk > initMus) initMuk = initMus;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getInitialMuStatic() const {
    return initMus;
}
//_____________________________________________________________________________
void
ExponentialSpringParameters::
setInitialMuKinetic(Real muk) {
    initMuk = muk;
    // Keep initMuk greater than or equal to 0.0.
    if(initMuk < 0.0) initMuk = 0.0;
    // Make sure initMus is greater than or equal to initMuk
    if(initMus < initMuk) initMus = initMuk;
}
//_____________________________________________________________________________
Real
ExponentialSpringParameters::
getInitialMuKinetic() const {
    return initMuk;
}

