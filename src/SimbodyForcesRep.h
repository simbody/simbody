#ifndef SimTK_SIMBODY_FORCES_REP_H_
#define SimTK_SIMBODY_FORCES_REP_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * Define the private implementations of some basic types of force
 * subsystems of use to Simbody users.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/MatterSubsystem.h"
#include "simbody/internal/SimbodyForces.h"

#include "ForceSubsystemRep.h"


namespace SimTK {

// 
// Define a linear spring between two stations s1 and s2 of a matter subsystem
// (a station is a point fixed on a particular body). A spring has a stiffness k, 
// and a natural length x0 at which it generates no force. Define the separation
// vector v=s2-s1, with x=|v| the spring's current length.
//
// We will request parameters in the State for k and x0 but require fixed stations.
// Defaults for k and x0 must be provided on construction.
//
// Then the potential energy stored in the spring is 
//    pe = k(x-x0)^2/2
// Forces are generated on both points, as the negative gradient of the
// potential energy at that point: 
//    f1 = d pe/d s1 =  k(x-x0)v/x
//    f2 = d pe/d s2 = -k(x-x0)v/x.
// Note that force is undefined when x=0; we'll return NaN vectors in that case.

class TwoPointSpringSubsystemRep : public ForceSubsystemRep {

    // state entries
    struct TwoPointLinearSpringParameters {
        TwoPointLinearSpringParameters()
          : body1(-1), body2(-1)
        { 
            station1.setToNaN(); station2.setToNaN();
            stiffness = naturalLength = CNT<Real>::getNaN();
        }

        TwoPointLinearSpringParameters(
            int b1, const Vec3& s1, int b2, const Vec3& s2,
            const Real& k, const Real& x0) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            stiffness(k), naturalLength(x0) 
        { 
            assert(b1 >= 0 && b2 >= 0 && b1 != b2);
            assert(stiffness >= 0 && naturalLength >= 0);
        }

        int  body1, body2;
        Vec3 station1, station2;    // in body frames
        Real stiffness, naturalLength;
    };

    struct GlobalMobilityDampingParameters {
        GlobalMobilityDampingParameters() : damping(0) { }
        GlobalMobilityDampingParameters(const Real& c)
          : damping(c) 
        { 
            assert(damping >= 0);
        }
        Real damping;   // 0 means "none"
    };

    struct Parameters {
        Parameters() : enabled(true) { }
        bool enabled;
        std::vector<TwoPointLinearSpringParameters>  twoPtLinearSprings;
        std::vector<GlobalMobilityDampingParameters> globalMobilityDampers;
    };

    // topological variables
    Parameters defaultParameters;

    // These must be filled in during realizeConstruction and treated
    // as const thereafter. These are garbage unless built=true.
    mutable int parameterVarsIndex;
    mutable bool built;

    const Parameters& getParameters(const State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,parameterVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,parameterVarsIndex)).upd();
    }

public:
    TwoPointSpringSubsystemRep()
     : ForceSubsystemRep("TwoPointSpringSubsystem", "0.0.1"), 
       built(false)
    {
    }

    int addLinearTwoPointSpring(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& stiffness,
                                const Real& naturalLength)
    {
        assert(stiffness >= 0);
        assert(naturalLength >= 0);
        assert(body1 != body2);
        defaultParameters.twoPtLinearSprings.push_back(
            TwoPointLinearSpringParameters(body1,s1,body2,s2,stiffness,naturalLength));
        return (int)defaultParameters.twoPtLinearSprings.size() - 1;
    }

    int addGlobalMobilityDamping(const Real& dampingFactor) {
        assert(dampingFactor >= 0);
        defaultParameters.globalMobilityDampers.push_back(
            GlobalMobilityDampingParameters(dampingFactor));
        return (int)defaultParameters.globalMobilityDampers.size() - 1;
    }

    void realizeConstruction(State& s) const {
        parameterVarsIndex = s.allocateDiscreteVariable(getMySubsystemIndex(), Stage::Parametrized, 
            new Value<Parameters>(defaultParameters));
        built = true;
    }

    void realizeModeling(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeParameters(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizeConfiguration(const State& s) const {
        // Nothing to compute here.
    }

    void realizeMotion(const State& s) const {
        // Nothing to compute here.
    }

    void realizeDynamics(const State& s) const {
        const Parameters& p = getParameters(s);
        if (!p.enabled) return;

        const MultibodySystem& mbs    = getMultibodySystem(); // my owner
        const MatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to system-global cache entries.
        Real&                  pe              = mbs.updPotentialEnergy(s);
        Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s);
        Vector&                mobilityForces  = mbs.updMobilityForces(s);

        // Global mobility dampers

        for (int i=0; i < (int)p.globalMobilityDampers.size(); ++i) {
            const Real c = p.globalMobilityDampers[i].damping;
            if (c == 0) continue;
            mobilityForces -= c*matter.getU(s);
        }

        // Linear two-point springs

        for (int i=0; i < (int)p.twoPtLinearSprings.size(); ++i) {
            const TwoPointLinearSpringParameters& spring =
                p.twoPtLinearSprings[i];
            const Transform& X_GB1 = matter.getBodyConfiguration(s, spring.body1);
            const Transform& X_GB2 = matter.getBodyConfiguration(s, spring.body2);

            const Vec3 s1_G = X_GB1.R() * spring.station1;
            const Vec3 s2_G = X_GB2.R() * spring.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            const Real d   = r_G.norm();  // distance between the points
            const Real stretch   = d - spring.naturalLength; // + -> tension, - -> compression
            const Real frcScalar = spring.stiffness*stretch; // k(x-x0)

            pe += 0.5 * frcScalar * stretch; // 1/2 k (x-x0)^2

            const Vec3 f1_G = (frcScalar/d) * r_G;
            rigidBodyForces[spring.body1] +=  SpatialVec(s1_G % f1_G, f1_G);
            rigidBodyForces[spring.body2] -=  SpatialVec(s2_G % f1_G, f1_G);
        }
    }

    void realizeReaction(const State& s) const {
        // Nothing to compute here.
    }

    TwoPointSpringSubsystemRep* cloneSubsystemRep() const {return new TwoPointSpringSubsystemRep(*this);}
    friend std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::Parameters&); 
};
// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::Parameters&) 
{assert(false);return o;}

// This is an empty placeholder force subsystem. It does nothing but exist; is that
// really so different from the rest of us?
class EmptyForcesSubsystemRep : public ForceSubsystemRep {
public:
    EmptyForcesSubsystemRep()
      : ForceSubsystemRep("EmptyForcesSubsystem", "0.0.1") { }

    EmptyForcesSubsystemRep* cloneSubsystemRep() const 
      { return new EmptyForcesSubsystemRep(*this); }

    SimTK_DOWNCAST(EmptyForcesSubsystemRep,ForceSubsystemRep);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_REP_H_
