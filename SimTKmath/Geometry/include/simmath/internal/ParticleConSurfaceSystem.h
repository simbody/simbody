#ifndef SimTK_SIMBODY_PARTICLECONSURFACESYSTEM_H_
#define SimTK_SIMBODY_PARTICLECONSURFACESYSTEM_H_

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
 * along a smooth surface. It is used to integrate Geodesic curves along
 * implicit surfaces. The particle on surface constraint is stabilized
 * with coordinate projection. (Implementation based on PendulumSystem.h)
 *
 **/

#include "simmath/internal/BicubicSurface.h" // XXX compiler needed this
#include "simmath/internal/ContactGeometry.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

namespace SimTK {

class ParticleConSurfaceSystem;
class ParticleConSurfaceSystemGuts: public System::Guts {
    friend class ParticleConSurfaceSystem;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
//    mutable DiscreteVariableIndex massIndex, lengthIndex, gravityIndex;
    DiscreteVariableIndex geodesicIndex, geometryIndex;
    mutable QIndex q0;
    mutable UIndex u0;
    mutable QErrIndex qerr0;
    mutable UErrIndex uerr0;
    mutable UDotErrIndex udoterr0;
    mutable EventTriggerByStageIndex event0;

public:
    ParticleConSurfaceSystemGuts(const ContactGeometryImpl& geom)
    : Guts(), geom(geom) {
        // Index types set themselves invalid on construction.
    }

    inline const ParticleConSurfaceSystem& getParticleConSurfaceSystem() const;
    
    SubsystemIndex getSubsysIndex() const {
        return subsysIndex;
    }

    /*virtual*/ParticleConSurfaceSystemGuts* cloneImpl() const override {return new ParticleConSurfaceSystemGuts(*this);}

        /////////////////////////////////////////////////////////
        // Implementation of continuous DynamicSystem virtuals //
        /////////////////////////////////////////////////////////

    /*virtual*/int realizeTopologyImpl(State&) const override;
    /*virtual*/int realizeModelImpl(State&) const override;
    /*virtual*/int realizeInstanceImpl(const State&) const override;
    /*virtual*/int realizePositionImpl(const State&) const override;
    /*virtual*/int realizeVelocityImpl(const State&) const override;
    /*virtual*/int realizeDynamicsImpl(const State&) const override;
    /*virtual*/int realizeAccelerationImpl(const State&) const override;

    // qdot==u here so these are just copies
    /*virtual*/void multiplyByNImpl(const State& state, const Vector& u, 
                                 Vector& dq) const override {dq=u;}
    /*virtual*/void multiplyByNTransposeImpl(const State& state, const Vector& fq, 
                                          Vector& fu) const override {fu=fq;}
    /*virtual*/void multiplyByNPInvImpl(const State& state, const Vector& dq, 
                                     Vector& u) const override {u=dq;}
    /*virtual*/void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu, 
                                              Vector& fq) const override {fq=fu;}

    // No prescribed motion.
    /*virtual*/bool prescribeQImpl(State&) const override {return false;}
    /*virtual*/bool prescribeUImpl(State&) const override {return false;}

    // No constraints.
    /*virtual*/void projectQImpl(State&, Vector& qErrEst,
             const ProjectOptions& options, ProjectResults& results) const override;
    /*virtual*/void projectUImpl(State&, Vector& uErrEst,
             const ProjectOptions& options, ProjectResults& results) const override;
private:
    const ContactGeometryImpl& geom;
}; // class ParticleConSurfaceSystemGuts



class ParticleConSurfaceSystem: public System {
public:
    ParticleConSurfaceSystem(const ContactGeometryImpl& geom) : System()
    { 
        adoptSystemGuts(new ParticleConSurfaceSystemGuts(geom));
        DefaultSystemSubsystem defsub(*this);
        updGuts().subsysIndex = defsub.getMySubsystemIndex();

        setHasTimeAdvancedEvents(false);
    }

    const ParticleConSurfaceSystemGuts& getGuts() const {
        return SimTK_DYNAMIC_CAST_DEBUG<const ParticleConSurfaceSystemGuts&>
                                                            (getSystemGuts());
    }

    ParticleConSurfaceSystemGuts& updGuts() {
        return SimTK_DYNAMIC_CAST_DEBUG<ParticleConSurfaceSystemGuts&>
                                                            (updSystemGuts());
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        const ParticleConSurfaceSystemGuts& guts = getGuts();
        updDefaultState().updU(guts.subsysIndex) = u;
        updDefaultState().updQ(guts.subsysIndex) = q;
        updDefaultState().updTime() = t;
    }


}; // class ParticleConSurfaceSystem


inline const ParticleConSurfaceSystem& ParticleConSurfaceSystemGuts::
getParticleConSurfaceSystem() const {
    return static_cast<const ParticleConSurfaceSystem&>(getSystem());
}


} // namespace SimTK

#endif /*SimTK_SIMBODY_PARTICLECONSURFACESYSTEM_H_*/
