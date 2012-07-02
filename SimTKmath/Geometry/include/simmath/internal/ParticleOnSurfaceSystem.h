#ifndef SimTK_SIMBODY_PARTICLEONGSURFACESYSTEM_H_
#define SimTK_SIMBODY_PARTICLEONGSURFACESYSTEM_H_

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

#include "simmath/internal/BicubicSurface.h" // XXX compiler needed this
#include "simmath/internal/ContactGeometry.h"
#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

namespace SimTK {

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
    ParticleOnSurfaceSystemGuts(const ContactGeometry& geom, double alpha, double beta)
    : Guts(), geom(geom), alpha(alpha), beta(beta) {
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
    Real alpha;
    Real beta;
}; // class ParticleOnSurfaceSystemGuts



class ParticleOnSurfaceSystem: public System {
public:
    ParticleOnSurfaceSystem(const ContactGeometry& geom, double alpha, double beta) : System()
    { 
        adoptSystemGuts(new ParticleOnSurfaceSystemGuts(geom, alpha, beta));
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




/*
 * This system is a 3d particle mass constrained to move along a surface
 * with no applied force (other than the constraint reaction force normal
 * to the surface). With no applied for the particle traces a geodesic along
 * the surface.
 *
 *
 * The DAE for a generic multibody system is:
 *       qdot = Qu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r be the 3d coordinates of the particle
 *       g(r) be the implicit surface function
 *       G be the gradient of the implicit surface function
 *       H be the hessian of the implicit surface function
 *
 * We will express implicit surface constraint as
 *                g(r) = 0    (perr)
 *                 Gr' = 0    (verr)
 *        Gr'' = -G'r' = - ~r'Hr' (aerr)
 *
 * So the matrix A = G and b = -r'Hr', and the
 * equations of motion are:
 *     [ M  ~G ] [ r'' ]   [     0    ]
 *     [ G   0 ] [  L  ] = [ - ~r'Hr' ]
 * where L (the Lagrange multiplier) is proportional to
 * the constraint force.
 *
 * solving for L,
 *        L = (GM\~G)\~r'Hr'
 *      r'' = - M\~G(GM\~G)\~r'Hr'
 *
 *      where ~A = transpose of A
 *        and A\ = inverse of A
 */
int ParticleOnSurfaceSystemGuts::realizeTopologyImpl(State& s) const {

    const Vector init(3, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeInstanceImpl(const State& s) const {

    System::Guts::realizeInstanceImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizePositionImpl(const State& s) const {

    System::Guts::realizePositionImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];
    qdot[2] = u[2];

    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeDynamicsImpl(const State& s) const {

    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeAccelerationImpl(const State& s) const {

    // XXX assume unit mass

    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       udot = s.updUDot(subsysIndex);

    Vec3 v(u[0], u[1], u[2]);
    Vec3 a(0);

    Real g = geom.calcSurfaceValue(q);
    Vec3 GT = geom.calcSurfaceGradient(q);
    Mat33 H = geom.calcSurfaceHessian(q);
    Real Gdotu = ~v*(H*v);
    Real L = (Gdotu + beta*~GT*v + alpha*g)/(~GT*GT);
    a = GT*-L;

    udot[0] = a[0]; udot[1] = a[1]; udot[2] = a[2];
    s.updQDotDot() = udot;

    System::Guts::realizeAccelerationImpl(s);
    return 0;
}


} // namespace SimTK

#endif /*SimTK_SIMBODY_PARTICLEONGSURFACESYSTEM_H_*/
