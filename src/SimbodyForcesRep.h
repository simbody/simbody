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
 * Define the private implementation of the SimbodyForcesSubsystem
 * class (a kind of MechanicalForcesSubsystem).
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
    struct Parameters {
        Parameters(const Real& k, const Real& x0) 
            : stiffness(k), naturalLength(x0), gravity(0), damping(0) { }
        Real stiffness, naturalLength;
        Vec3 gravity;
        Real damping;
    };

    struct ConfigurationCache {
        Vec3 station1_G, station2_G; // body station vectors re-expressed in G
        Vec3 v_G;
        Real x;       // length of above vector
        Real fscalar; // k(x-x0)
        Real pe;
    };
    struct DynamicsCache {
        Vec3 f1_G;    // f2 is the negative of this
    };

    // topological variables
    int  body1, body2;
    Vec3 station1, station2;
    Parameters defaultParameters;

    // These must be filled in during realizeConstruction and treated
    // as const thereafter. These are garbage unless built=true.
    mutable int parameterVarsIndex;
    mutable int configurationCacheIndex;
    mutable int dynamicsCacheIndex;
    mutable bool built;

    const Parameters& getParameters(const State& s) const {
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,parameterVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,parameterVarsIndex)).upd();
    }
    const ConfigurationCache& getConfigurationCache(const State& s) const {
        return Value<ConfigurationCache>::downcast(
            getCacheEntry(s,configurationCacheIndex)).get();
    }
    ConfigurationCache& updConfigurationCache(const State& s) const {
        return Value<ConfigurationCache>::downcast(
            updCacheEntry(s,configurationCacheIndex)).upd();
    }
    const DynamicsCache& getDynamicsCache(const State& s) const {
        return Value<DynamicsCache>::downcast(
            getCacheEntry(s,dynamicsCacheIndex)).get();
    }
    DynamicsCache& updDynamicsCache(const State& s) const {
        return Value<DynamicsCache>::downcast(
            updCacheEntry(s,dynamicsCacheIndex)).upd();
    }


public:
    TwoPointSpringSubsystemRep(int b1, const Vec3& s1,
                               int b2, const Vec3& s2,
                               const Real& k, const Real& x0)
     : ForceSubsystemRep("TwoPointSpringSubsystem", "0.0.1"), 
       body1(b1), body2(b2), station1(s1), station2(s2),
       defaultParameters(k,x0), built(false)
    {
    }

    // Pure virtuals
    /// This is a Configured stage operator.
    Real calcPotentialEnergy(const State& s) const {
        SimTK_STAGECHECK_GE(getStage(s), Stage::Configured, 
            "TwoPointSpringSubsystem::calcPotentialEnergy()");
        return getConfigurationCache(s).pe;
    }

    /// This is a Dynamics stage operator.
    void addInForces(const State& s, const MatterSubsystem& matter,
                     Vector_<SpatialVec>& rigidBodyForces,
                     Vector_<Vec3>&       particleForces,
                     Vector&              mobilityForces) const 
    {
        assert(matter.getMySubsystemIndex() == getMatterSubsystemIndex());
        SimTK_STAGECHECK_GE(getStage(s), Stage::Dynamics, 
            "TwoPointSpringSubsystem::addInForces()");

        const ConfigurationCache& cc = getConfigurationCache(s);
        const DynamicsCache&      dc = getDynamicsCache(s);
        rigidBodyForces[body1] += SpatialVec( cc.station1_G % dc.f1_G, dc.f1_G );
        rigidBodyForces[body2] -= SpatialVec( cc.station2_G % dc.f1_G, dc.f1_G );

        if (getGravity(s) != Vec3(0))
            matter.addInGravity(s, getGravity(s), rigidBodyForces);

        if (getDamping(s) != 0)
            mobilityForces += -getDamping(s)*matter.getU(s);
    }

    const Vec3& getGravity(const State& s) const {return getParameters(s).gravity;}
    Vec3&       updGravity(State& s)       const {return updParameters(s).gravity;}

    const Real& getDamping(const State& s) const {return getParameters(s).damping;}
    Real&       updDamping(State& s)       const {return updParameters(s).damping;}

    const Real& getStiffness(const State& s) const {return getParameters(s).stiffness;}
    Real&       updStiffness(State& s)       const {return updParameters(s).stiffness;}

    const Real& getNaturalLength(const State& s) const {return getParameters(s).naturalLength;}
    Real&       updNaturalLength(State& s)       const {return updParameters(s).naturalLength;}

    const Real& getPotentialEnergy(const State& s) const {return getConfigurationCache(s).pe;}
    const Vec3& getForceOnStation1(const State& s) const {return getDynamicsCache(s).f1_G;}

    void realizeConstruction(State& s) const {
        parameterVarsIndex = s.allocateDiscreteVariable(getMySubsystemIndex(), Stage::Parametrized, 
            new Value<Parameters>(defaultParameters));
        configurationCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), Stage::Configured,
            new Value<ConfigurationCache>());
        dynamicsCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), Stage::Dynamics,
            new Value<DynamicsCache>());
        built = true;
    }

    void realizeModeling(State& s) const {
        static const char* loc = "TwoPointSpring::realizeModeling()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Built, loc);
        // Sorry, no choices available at the moment.
    }

    void realizeParameters(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeParameters()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Modeled, loc);
        SimTK_VALUECHECK_NONNEG_ALWAYS(getStiffness(s), "stiffness", loc);
        SimTK_VALUECHECK_NONNEG_ALWAYS(getNaturalLength(s), "naturalLength", loc);

        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeTime()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Parametrized, loc);
        // Nothing to compute here.
    }

    void realizeConfiguration(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeConfiguration()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Timed, loc);

        const Parameters&   p  = getParameters(s);
        ConfigurationCache& cc = updConfigurationCache(s);

        const Transform& X_GB1 = getMatterSubsystem().getBodyConfiguration(s, body1);
        const Transform& X_GB2 = getMatterSubsystem().getBodyConfiguration(s, body2);

        // Fill in the configuration cache.
        cc.station1_G = X_GB1.R() * station1;   // stations expressed in G will be needed later
        cc.station2_G = X_GB2.R() * station2;

        const Vec3 p1_G = X_GB1.T() + cc.station1_G;    // station point locations in ground
        const Vec3 p2_G = X_GB2.T() + cc.station2_G;

        cc.v_G              = p2_G - p1_G;
        cc.x                = cc.v_G.norm();
        const Real stretch  = cc.x - p.naturalLength;   // + -> tension, - -> compression
        cc.fscalar          = p.stiffness * stretch;    // k(x-x0)
        cc.pe               = 0.5 * cc.fscalar * stretch;
    }

    void realizeMotion(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeMotion()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Configured, loc);
        // Nothing to compute here.
    }

    void realizeDynamics(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeDynamics()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Moving, loc);

        const ConfigurationCache& cc = getConfigurationCache(s);
        DynamicsCache&            dc = updDynamicsCache(s);

        dc.f1_G = (cc.fscalar/cc.x) * cc.v_G;  // NaNs if x (and hence v) is 0
    }

    void realizeReaction(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeReaction()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Dynamics, loc);

        // Nothing to compute here.
    }

    TwoPointSpringSubsystemRep* cloneSubsystemRep() const {return new TwoPointSpringSubsystemRep(*this);}
    friend std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::Parameters&); 
    friend std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::ConfigurationCache&); 
    friend std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::DynamicsCache&);

};
// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::Parameters&) 
{assert(false);return o;}
std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::ConfigurationCache&) 
{assert(false);return o;}
std::ostream& operator<<(std::ostream& o, 
                         const TwoPointSpringSubsystemRep::DynamicsCache&) 
{assert(false);return o;}

class EmptyForcesSubsystemRep : public ForceSubsystemRep {
public:
    EmptyForcesSubsystemRep()
      : ForceSubsystemRep("EmptyForcesSubsystem", "0.0.1") { }
    EmptyForcesSubsystemRep(const MatterSubsystem& m) 
      : ForceSubsystemRep("EmptyForcesSubsystem", "0.0.1") { 
        setMatterSubsystemIndex(m.getMySubsystemIndex());
    }


    // Pure virtuals
    // This is a Configured stage operator.
    Real calcPotentialEnergy(const State& s) const {
        SimTK_STAGECHECK_GE(getStage(s), Stage::Configured, 
            "EmptyForcesSubsystem::calcPotentialEnergy()");
        return 0;
    }

    // This is a Dynamics stage operator.
    void addInForces(const State& s, const MatterSubsystem& matter,
                     Vector_<SpatialVec>& rigidBodyForces,
                     Vector_<Vec3>&       particleForces,
                     Vector&              mobilityForces) const 
    {
        assert(matter.getMySubsystemIndex() == getMatterSubsystemIndex());
        SimTK_STAGECHECK_GE(getStage(s), Stage::Dynamics, 
            "EmptyForcesSubsystem::addInForces()");
        // nothing to add
    }

    EmptyForcesSubsystemRep* cloneSubsystemRep() const 
      { return new EmptyForcesSubsystemRep(*this); }

    SimTK_DOWNCAST(EmptyForcesSubsystemRep,ForceSubsystemRep);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_REP_H_
