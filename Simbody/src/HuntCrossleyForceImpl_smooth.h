#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_SMOOTH_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_SMOOTH_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Antoine Falisse, Gil Serrancoli                              *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/HuntCrossleyForce_smooth.h"
#include "ForceImpl.h"

namespace SimTK {

class HuntCrossleyForceImpl_smooth : public ForceImpl {
public:
    class Parameters {
    public:
        Parameters() : stiffness(0), dissipation(0), staticFriction(0),
            dynamicFriction(0), viscousFriction(0), transitionVelocity(0) {
        }
        Parameters(Real stiffness, Real dissipation, Real staticFriction,
            Real dynamicFriction, Real viscousFriction,
            Real transitionVelocity) : stiffness(stiffness),
            dissipation(dissipation), staticFriction(staticFriction),
            dynamicFriction(dynamicFriction), viscousFriction(viscousFriction),
            transitionVelocity(transitionVelocity) {
        }
        Real stiffness, dissipation, staticFriction, dynamicFriction,
            viscousFriction, transitionVelocity;
    };
    struct InstanceVars {
        InstanceVars(const MobilizedBody defbody, const Real defradius,
            const Vec3 defLocContactSphere, const Parameters defparameters)
            : BodySphere(defbody), RadiusContactSphere(defradius),
            LocContactSphere(defLocContactSphere), parameters(defparameters) {
        }
        MobilizedBody BodySphere;
        Real RadiusContactSphere;
        Vec3 LocContactSphere;
        Parameters parameters;
    };

    Real            RadiusContactSphere;
    Vec3            LocContactSphere;
    MobilizedBody   BodySphere;
    Parameters      parameters;
    Plane           ContactPlane;

    HuntCrossleyForceImpl_smooth(GeneralForceSubsystem& subsystem);

    HuntCrossleyForceImpl_smooth* clone() const override {
        return new HuntCrossleyForceImpl_smooth(*this);
    }
    /**
     * Set the contact material parameters
     *
     * @param stiffness       the stiffness constant
     * @param dissipation     the dissipation coefficient (c)
     * @param staticFriction  the coefficient of static friction (us)
     * @param dynamicFriction the coefficient of dynamic friction (ud)
     * @param viscousFriction the coefficient of viscous friction (uv)
     */
    void setParameters(Real stiffness, Real dissipation, Real staticFriction,
        Real dynamicFriction, Real viscousFriction, Real transitionVelocity);
    /** Get parameters */
    const Parameters& getParameters() const;
    /** Update parameters */
    Parameters& updParameters();
    /** Set the stiffness constant */
    void setStiffness(Real stiffness);
    /** Set the dissipation coefficient */
    void setDissipation(Real dissipation);
    /** Set the coefficient of static friction */
    void setStaticFriction(Real staticFriction);
    /** Set the coefficient of dynamic friction */
    void setDynamicFriction(Real dynamicFriction);
    /** Set the coefficient of viscous friction */
    void setViscousFriction(Real viscousFriction);
    /** Set the transition velocity */
    void setTransitionVelocity(Real transitionVelocity);
    /**
    * Set contact plane
    *
    * @param normal     direction of the normal to the plane of contact
    * @param offset     distance to the ground origin along the normal
    */
    void setContactPlane(Vec3 normal, Real offset);
    /** Set the Mobilized Body to which the contact sphere is attached */
    void setContactSphere(MobilizedBody bodyInput);
    /** Set the location of the contact sphere in the body frame */
    void setLocContactSphere(Vec3 LocContactSphere);
    /** Set the radius of the contact sphere */
    void setRadiusContactSphere(Real radius);
    /** Get the Mobilized Body to which the contact sphere is attached */
    MobilizedBody getBodySphere();
    /** Get the location of the contact sphere in the body frame */
    Vec3 getLocContactSphere();
    /** Set the radius of the sphere */
    Real getRadiusContactSphere();
    /** Get the location of the contact point in the body frame */
    Vec3 getContactPointInBody(const State& state);
    /** Get the location of the contact point in the ground frane */
    void getContactPointSphere(const State& state,Vec3& contactPointPos) const;
    /** Calculate contact force */
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>& particleForces, Vector& mobilityForces) const override;
    /** TODO */
    Real calcPotentialEnergy(const State& state) const override;
    void realizeTopology(State& state) const override;
private:
    const GeneralForceSubsystem&          subsystem;
    mutable CacheEntryIndex               energyCacheIndex;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_SMOOTH_H_
