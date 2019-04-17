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
            const Vec3 deflocSphere, const Parameters defparameters)
			: BodySphere(defbody), RadiusSphere(defradius),
            LocSphere(deflocSphere), parameters(defparameters) {
        }
		MobilizedBody BodySphere;
		Real RadiusSphere;
		Vec3 LocSphere;
		Parameters parameters;
	};

	Real            RadiusSphere;
	Vec3            LocSphere;
	MobilizedBody   BodySphere;
	Parameters		parameters;
	Plane			GroundPlane;

	HuntCrossleyForceImpl_smooth(GeneralForceSubsystem& subsystem);

	HuntCrossleyForceImpl_smooth* clone() const {
		return new HuntCrossleyForceImpl_smooth(*this);
	}

    /** Set parameters */
    void setParameters(Real stiffness, Real dissipation, Real staticFriction,
        Real dynamicFriction, Real viscousFriction, Real transitionVelocity);
    /** Get parameters */
    const Parameters& getParameters() const;
    /** Update parameters */
    Parameters& updParameters();
    /** Set stiffness */
    void setStiffness(Real stiffness);
    /** Set dissipation */
    void setDissipation(Real dissipation);
    /** Set static friction */
    void setStaticFriction(Real staticFriction);
    /** Set dynamic friction */
    void setDynamicFriction(Real dynamicFriction);
    /** Set viscous friction */
    void setViscousFriction(Real viscousFriction);
    /** Set transition velocity */
    void setTransitionVelocity(Real transitionVelocity);
	/** Set ground plane */
	void setGroundPlane(Vec3 normal, Real offset);
    /** Set Mobilized Body the sphere is attached to */
	void setBodySphere(MobilizedBody bodyInput);
	/** Set location of the sphere in the body frame */
	void setLocSphere(Vec3 locSphere);
	/** Set radius of the sphere */
	void setRadiusSphere(Real radius);
	/** Get Mobilized Body the sphere is attached to */
    MobilizedBody getBodySphere();
	/** Get location of the sphere in the body frame*/
	Vec3 getLocSphere();
	/** Get radius of the sphere */
	Real getRadiusSphere();
    /** Get location of the contact point in the body frame */
	Vec3 getContactPointInBody(const State& state);
    /** Get location of the contact point in the ground frane */
    void getContactPointSphere(const State& state,Vec3& contactPointPos) const;
	/** Calculate contact force */
	void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
		Vector_<Vec3>& particleForces, Vector& mobilityForces) const override;
    /** TODO */
	Real calcPotentialEnergy(const State& state) const;
    void realizeTopology(State& state) const;
private:
	const GeneralForceSubsystem&          subsystem;
	mutable CacheEntryIndex               energyCacheIndex;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_SMOOTH_H_
