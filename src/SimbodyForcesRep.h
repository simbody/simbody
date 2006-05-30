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
#include "simbody/internal/SimbodyForces.h"

#include "MultibodySystemRep.h"


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

class TwoPointSpringSubsystemRep : public MechanicalForcesSubsystemRep {

    // state entries
    struct Parameters {
        Parameters(const Real& k, const Real& x0) : stiffness(k), naturalLength(x0) { }
        Real stiffness, naturalLength;
    };

    struct ConfigurationCache {
        Vec3 v_G;
        Real x;    // length of above vector
        Real fmag; // k(x-x0)
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

    const Stage& getStage(const State& s) const {
        return s.getStage(getMySubsystemIndex());
    }
    void invalidateStage(const State& s, Stage g) const {
        s.invalidateStage(getMySubsystemIndex(),g);
    }
    void advanceToStage(const State& s, Stage g) const {
        s.advanceToStage(getMySubsystemIndex(),g);
    }

    const Parameters& getParameters(const State& s) const {
        return Value<Parameters>::downcast(
            s.getDiscreteVariable(getMySubsystemIndex(),parameterVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        return Value<Parameters>::downcast(
            s.updDiscreteVariable(getMySubsystemIndex(),parameterVarsIndex)).upd();
    }
    const ConfigurationCache& getConfigurationCache(const State& s) const {
        return Value<ConfigurationCache>::downcast(
            s.getCacheEntry(getMySubsystemIndex(),configurationCacheIndex)).get();
    }
    ConfigurationCache& updConfigurationCache(const State& s) const {
        return Value<ConfigurationCache>::downcast(
            s.updCacheEntry(getMySubsystemIndex(),configurationCacheIndex)).upd();
    }
    const DynamicsCache& getDynamicsCache(const State& s) const {
        return Value<DynamicsCache>::downcast(
            s.getCacheEntry(getMySubsystemIndex(),dynamicsCacheIndex)).get();
    }
    DynamicsCache& updDynamicsCache(const State& s) const {
        return Value<DynamicsCache>::downcast(
            s.updCacheEntry(getMySubsystemIndex(),dynamicsCacheIndex)).upd();
    }


public:
    TwoPointSpringSubsystemRep(const MechanicalSubsystem& m,
                   int b1, const Vec3& s1,
                   int b2, const Vec3& s2,
                   const Real& k, const Real& x0)
     : MechanicalForcesSubsystemRep("TwoPointSpringSubsystem", "0.0.1", m), 
       body1(b1), body2(b2), station1(s1), station2(s2),
       defaultParameters(k,x0), built(false)
    {
    }

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
        advanceToStage(s, Stage::Built);
    }

    void realizeModeling(State& s) const {
        static const char* loc = "TwoPointSpring::realizeModeling()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Built, loc);
        // Sorry, no choices available at the moment.
        advanceToStage(s, Stage::Modeled);
    }

    void realizeParameters(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeParameters()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Modeled, loc);
        SimTK_VALUECHECK_NONNEG_ALWAYS(getStiffness(s), "stiffness", loc);
        SimTK_VALUECHECK_NONNEG_ALWAYS(getNaturalLength(s), "naturalLength", loc);

        // Nothing to compute here.
        advanceToStage(s, Stage::Parametrized);
    }

    void realizeTime(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeTime()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Parametrized, loc);
        // Nothing to compute here.
        advanceToStage(s, Stage::Timed);
    }

    void realizeConfiguration(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeConfiguration()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Timed, loc);

        const Parameters&   p  = getParameters(s);
        ConfigurationCache& cc = updConfigurationCache(s);

        const Transform& X_GB1 = getMechanicalSubsystem().getBodyConfiguration(s, body1);
        const Transform& X_GB2 = getMechanicalSubsystem().getBodyConfiguration(s, body2);
        const Vec3 p1_G = X_GB1*station1;
        const Vec3 p2_G = X_GB2*station2;

        // Fill in the configuration cache.
        cc.v_G              = p2_G - p1_G;
        cc.x                = cc.v_G.norm();
        const Real stretch  = cc.x - p.naturalLength;   // + -> tension, - -> compression
        cc.fmag             = p.stiffness * stretch;    // k(x-x0)
        cc.pe               = 0.5 * cc.fmag * stretch;

        advanceToStage(s, Stage::Configured);
    }

    void realizeMotion(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeMotion()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Configured, loc);
        // Nothing to compute here.
        advanceToStage(s, Stage::Moving);
    }

    void realizeDynamics(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeDynamics()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Moving, loc);

        const ConfigurationCache& cc = getConfigurationCache(s);
        DynamicsCache&            dc = updDynamicsCache(s);

        dc.f1_G = (cc.fmag/cc.x) * cc.v_G;  // NaNs if x (and hence v) is 0

        advanceToStage(s, Stage::Dynamics);
    }

    void realizeReaction(const State& s) const {
        static const char* loc = "TwoPointSpring::realizeReaction()";
        SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Dynamics, loc);
        // Nothing to compute here.
        advanceToStage(s, Stage::Reacting);
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

/*
class SimbodyForcesSubsystemRep : public MechanicalForcesSubsystemRep {
    struct ForcesCache {
        Vector              mobilityForces;
        Vector_<SpatialVec> spatialForces;
    };
    const ForcesCache& getForcesCache(const State& s) const {
        return Value<ForcesCache>::downcast
            (s.getCacheEntry(forcesCacheIndex)).get();
    }
    ForcesCache& updForcesCache(const State& s) const {
        return Value<ForcesCache>::downcast
            (s.updCacheEntry(forcesCacheIndex)).upd();
    }
    const Real& getPotentialEnergyCache(const State& s) const {
        return Value<Real>::downcast
            (s.getCacheEntry(potentialEnergyCacheIndex)).get();
    }
    Real& updPotentialEnergyCache(const State& s) const {
        return Value<Real>::downcast
            (s.updCacheEntry(potentialEnergyCacheIndex)).upd();
    } 

public:
    SimbodyForcesSubsystemRep(const MechanicalSubsystem& m) //TODO
      : MechanicalForcesSubsystemRep(m) 
    { 
    }

    // This response is available at Stage::Configured.
    const Real& getPotentialEnergy(const State& s) const {
        return getPotentialEnergyCache(s);
    }

    // These responses are available at Stage::Moving
    const Vector& getMobilityForces(const State& s) const {
        const ForcesCache& f = getForcesCache(s);
        return f.mobilityForces;
    }
    const Vector_<SpatialVec>& getSpatialForces(const State& s) const {
        const ForcesCache& f = getForcesCache(s);
        return f.spatialForces;
    }

    // Construction: add force elements.
    void addGravity(const Vec3& g);
    void addMobilitySpring(int body, int dof,
                           const Real& stiffness,
                           const Real& naturalLength);
    void addTwoPointSpring(int body1, const Vec3& s1,
                           int body2, const Vec3& s2,
                           const Real&    stiffness,
                           const Real&    naturalLength);
    void addTwoPointDamper(int body1, const Vec3& s1,
                           int body2, const Vec3& s2,
                           const Real&    damping);

    int getNElements() const {return elements.size();}

    SimbodyForcesSubsystemRep* cloneSubsystemRep() const 
      { return new SimbodyForcesSubsystemRep(*this); }

    void realizeConstruction(State& s) const {
        // These indices are mutable here but must be treated as const everywhere else.
        potentialEnergyCacheIndex = 
            s.allocateCacheEntry(Stage::Configured, new Value<Real>(CNT<Real>::getNaN()));
        forcesCacheIndex = 
            s.allocateCacheEntry(Stage::Moving, new Value<ForcesCache>());
    }
    void realizeModeling(State& s) const { 
        ForcesCache& f = updForcesCache(s);
        f.mobilityForces.resize(getMechanicalSubsystem().getNMobilities());
        f.mobilityForces.setToNaN();
        f.spatialForces.resize(getMechanicalSubsystem().getNBodies());
        f.spatialForces.setToNaN();
    }

    void realizeParameters   (const State&) const { }
    void realizeTime         (const State&) const { }
    void realizeConfiguration(const State& s) const {
        // Time to calculate potential energy
        Real pe = 0;
        for (int i=0; i < getNElements(); ++i) {
            const SimbodyForceElement& e = elements[i];
            if (!e.isPotentialForce())
                continue;
            pe += e.getPotentialEnergy(s);
        }
        updPotentialEnergyCache(s) = pe; 
    }
    void realizeMotion       (const State&) const { }
    void realizeDynamics     (const State&) const { 
    }
    void realizeReaction     (const State&) const { }

    SimTK_DOWNCAST(SimbodyForcesSubsystemRep,MechanicalForcesSubsystemRep);
private:
    // This is like the state cache for use during construction. These
    // must be treated as const after realizeConstruction().
    mutable int forcesCacheIndex;
    mutable int potentialEnergyCacheIndex;

    List<SimbodyForceElement> elements;
};
*/


class EmptyForcesSubsystemRep : public MechanicalForcesSubsystemRep {
public:
    EmptyForcesSubsystemRep()
      : MechanicalForcesSubsystemRep("EmptyForcesSubsystem", "0.0.1", MechanicalSubsystem()) { }
    EmptyForcesSubsystemRep(const MechanicalSubsystem& m) 
      : MechanicalForcesSubsystemRep("EmptyForcesSubsystem", "0.0.1", m) { }

    EmptyForcesSubsystemRep* cloneSubsystemRep() const 
      { return new EmptyForcesSubsystemRep(*this); }

    void realizeConstruction(State&) const { }
    void realizeModeling    (State&) const { }

    SimTK_DOWNCAST(EmptyForcesSubsystemRep,MechanicalForcesSubsystemRep);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_REP_H_
