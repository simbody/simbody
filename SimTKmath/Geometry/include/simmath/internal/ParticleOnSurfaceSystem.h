#ifndef SimTK_SIMBODY_PARTICLEONSURFACESYSTEM_H_
#define SimTK_SIMBODY_PARTICLEONSURFACESYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKmath                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness, Michael Sherman                                     *
 * Contributors:                                                              *
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


/**
 * This is a System that represents the dynamics of a particle moving
 * along a smooth surfaces. It is used to integrate Geodesic curves along
 * implicit surfaces. The particle on surface constraint is stabilized
 * with Baumgarte stabilization. (Implementation based on Pendulum system.)
 *
 **/

#include "SimTKcommon.h"
#include "simmath/internal/BicubicSurface.h" // XXX compiler needed this
#include "simmath/internal/ContactGeometry.h"

namespace SimTK {

const Real gamma = 0.01; // Baumgarte stabilization; params from Ascher1993
const Real alpha = gamma*gamma; // position stabilization
const Real beta = 2*gamma; // velocity stabilization

class ParticleOnSurfaceSystem;
class ParticleOnSurfaceSystemGuts: public System::Guts {
    friend class ParticleOnSurfaceSystem;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
//    mutable DiscreteVariableIndex massIndex, lengthIndex, gravityIndex;
    DiscreteVariableIndex geodesicIndex, geometryIndex;
    mutable QIndex q0;
    mutable UIndex u0;
//    mutable QErrIndex qerr0;
//    mutable UErrIndex uerr0;
//    mutable UDotErrIndex udoterr0;
//    mutable EventTriggerByStageIndex event0;
//    mutable CacheEntryIndex mgForceIndex; // a cache entry m*g calculated at Dynamics stage
public:
    ParticleOnSurfaceSystemGuts(const ContactGeometry& geom)
    : Guts(), geom(geom) {
        // Index types set themselves invalid on construction.
    }

    const ParticleOnSurfaceSystem& getParticleOnSurfaceSystem() const {
        return reinterpret_cast<const ParticleOnSurfaceSystem&>(getSystem());
    }
    
    SubsystemIndex getSubsysIndex() const {
        return subsysIndex;
    }

    /*virtual*/ParticleOnSurfaceSystemGuts* cloneImpl() const {return new ParticleOnSurfaceSystemGuts(*this);}

        /////////////////////////////////////////////////////////
        // Implementation of continuous DynamicSystem virtuals //
        /////////////////////////////////////////////////////////

    /*virtual*/int realizeTopologyImpl(State&) const;
    /*virtual*/int realizeModelImpl(State&) const;
    /*virtual*/int realizeInstanceImpl(const State&) const;
    /*virtual*/int realizePositionImpl(const State&) const;
    /*virtual*/int realizeVelocityImpl(const State&) const;
    /*virtual*/int realizeDynamicsImpl(const State&) const;
    /*virtual*/int realizeAccelerationImpl(const State&) const;

    // qdot==u here so these are just copies
    /*virtual*/void multiplyByNImpl(const State& state, const Vector& u, 
                                 Vector& dq) const {dq=u;}
    /*virtual*/void multiplyByNTransposeImpl(const State& state, const Vector& fq, 
                                          Vector& fu) const {fu=fq;}
    /*virtual*/void multiplyByNPInvImpl(const State& state, const Vector& dq, 
                                     Vector& u) const {u=dq;}
    /*virtual*/void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu, 
                                              Vector& fq) const {fq=fu;}

    // No prescribed motion.
//    /*virtual*/bool prescribeQImpl(State&) const {return false;}
//    /*virtual*/bool prescribeUImpl(State&) const {return false;}

    // No constraints.
//    /*virtual*/void projectQImpl(State&, Vector& qErrEst,
//             const ProjectOptions& options, ProjectResults& results) const {return;}
//    /*virtual*/void projectUImpl(State&, Vector& uErrEst,
//             const ProjectOptions& options, ProjectResults& results) const {return;}
private:
    ContactGeometry geom;

}; // class ParticleOnSurfaceSystemGuts



class ParticleOnSurfaceSystem: public System {
public:
    ParticleOnSurfaceSystem(const ContactGeometry& geom) : System()
    { 
        adoptSystemGuts(new ParticleOnSurfaceSystemGuts(geom));
        DefaultSystemSubsystem defsub(*this);
        updGuts().subsysIndex = defsub.getMySubsystemIndex();

//        setHasTimeAdvancedEvents(false);
    }

    const ParticleOnSurfaceSystemGuts& getGuts() const {
        return dynamic_cast<const ParticleOnSurfaceSystemGuts&>(getSystemGuts());
    }

    ParticleOnSurfaceSystemGuts& updGuts() {
        return dynamic_cast<ParticleOnSurfaceSystemGuts&>(updSystemGuts());
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        const ParticleOnSurfaceSystemGuts& guts = getGuts();
        updDefaultState().updU(guts.subsysIndex) = u;
        updDefaultState().updQ(guts.subsysIndex) = q;
        updDefaultState().updTime() = t;
    }


}; // class ParticleOnSurfaceSystem



} // namespace SimTK

#endif /*SimTK_SIMBODY_PARTICLEONSURFACESYSTEM_H_*/
