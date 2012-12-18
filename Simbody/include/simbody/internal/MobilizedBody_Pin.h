#ifndef SimTK_SIMBODY_MOBILIZED_BODY_PIN_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_PIN_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy, Peter Eastman                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
Declares the MobilizedBody::Pin class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// One mobility -- rotation about the common z axis of the inboard
/// and outboard mobilizer frames.
/// Synonym: Torsion
class SimTK_SIMBODY_EXPORT MobilizedBody::Pin : public MobilizedBody {
public:
        // SPECIALIZED INTERFACE FOR PIN MOBILIZER

    // "Angle" is just a nicer name for a pin joint's lone generalized coordinate q.
    Pin& setDefaultAngle(Real angleInRadians) {return setDefaultQ(angleInRadians);}
    Real getDefaultAngle() const              {return getDefaultQ();}

        // Friendly, mobilizer-specific access to generalized coordinates and speeds.

    void setAngle(State& s, Real angleInRadians) {setQ(s, angleInRadians);}
    Real getAngle(const State& s) const {return getQ(s);}

    void setRate(State& s, Real rateInRadiansPerTime) {setU(s, rateInRadiansPerTime);}
    Real getRate(const State& s) const {return getU(s);}

    // Mobility forces are "u-like", that is, one per dof.
    Real getAppliedPinTorque(const State& s, const Vector& mobilityForces) const {
        return getMyPartU(s,mobilityForces);
    }
    void applyPinTorque(const State& s, Real torque, Vector& mobilityForces) const {
        updMyPartU(s,mobilityForces) += torque;
    }

        // STANDARDIZED MOBILIZED BODY INTERFACE

        // required constructors
    explicit Pin(Direction=Forward);
    Pin(MobilizedBody& parent, const Body&, Direction=Forward);
    Pin(MobilizedBody& parent, const Transform& inbFrame,
        const Body&,           const Transform& outbFrame,
        Direction=Forward);

        // access to generalized coordinates q and generalized speeds u
    Pin& setDefaultQ(Real);
    Real getDefaultQ() const;

    Real getQ(const State&) const;
    Real getQDot(const State&) const;
    Real getQDotDot(const State&) const;
    Real getU(const State&) const;
    Real getUDot(const State&) const;

    void setQ(State&, Real) const;
    void setU(State&, Real) const;

    Real getMyPartQ(const State&, const Vector& qlike) const;
    Real getMyPartU(const State&, const Vector& ulike) const;
   
    Real& updMyPartQ(const State&, Vector& qlike) const;
    Real& updMyPartU(const State&, Vector& ulike) const;

        // specialize return type for convenience
    Pin& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this; }
    Pin& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g)
      { (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this; }
    Pin& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this; }
    Pin& setDefaultInboardFrame(const Transform& X_PF)
      { (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this; }
    Pin& setDefaultOutboardFrame(const Transform& X_BM)
      { (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this; }

    /** @cond **/ // Don't let doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Pin, PinImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_PIN_H_



