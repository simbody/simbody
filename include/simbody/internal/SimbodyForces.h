#ifndef SimTK_SIMBODY_FORCES_H_
#define SimTK_SIMBODY_FORCES_H_

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"

#include "simbody/internal/SimbodyState.h"

#include <cassert>

namespace SimTK {

class BasicMechanicalForceElements : public MechanicalForcesSubsystem {
public:
    BasicMechanicalForceElements(const MechanicalSubsystem& mech) 
        : MechanicalForcesSubsystem(mech), defaultGravity(0)
    {
    }

    MechanicalForcesSubsystem* cloneMechanicalForcesSubsystem() const {
        return new BasicMechanicalForceElements(*this);
    }

    void realizeConstruction(State&) { }
    void realizeModeling(State&) const { }

    void setGravity(const Vec3& g) { defaultGravity=g; }

    SimTK_DOWNCAST2(BasicMechanicalForceElements, MechanicalForcesSubsystem, Subsystem);
private:
    Vec3 defaultGravity;
};


// Construction invariants are:
//    body 1, body 2
//    default values for stations, naturalLength, stiffness
//    parametrized or not, discrete variable index if parameterized
//    cache index for results
// Modeling options are: 
//    none
// Parameters are:
//    none if (not parametrized)
//    otherwise: station locations, naturalLength, stiffness
// 
class LinearSpring {
    class Parameters;
    class Results;

public:
    LinearSpring(const SBStation& s1, const SBStation& s2, const Real& length0, const Real& k,
                 bool isParametrized=false)
      : body1(s1.body), body2(s2.body), 
        defaultState(s1.station, s2.station, length0, k), 
        parametrized(isParametrized), stateIndex(-1), cacheIndex(-1)
    {
        SimTK_ASSERT(body1 != body2, "LinearSpring()");
        SimTK_ASSERT(length0 >= 0., "LinearSpring()");
        SimTK_ASSERT(k >= 0., "LinearSpring()");
    }

    void realizeConstruction(State& s) {
        if (parametrized)
            stateIndex = s.allocateDiscreteVariable
                        (Stage::Parametrized,new Value<Parameters>(defaultState));
        cacheIndex = s.allocateCacheEntry
                     (Stage::Configured, new Value<Results>());
    }

    void realizeModeling(State& s) const { }

    void realizeParameters(const State& s) const {
        // TODO: if (parametrized)
        // TODO:   check parameter validity
        // TODO:   throw exception if bad
    }

    void realizeTime         (const State&) const { }

    void realizeConfiguration(const State& s, const SimbodySubsystem& t) const {
        Results&          c = updCacheEntries(s);
        const Parameters& v = parametrized ? getStateVariables(s) : defaultState;

        const Vec3 from1to2 = 
            t.getStationLocation(s, body2, v.station2) - t.getStationLocation(s, body1, v.station1);
        c.distance = from1to2.norm();
        c.fromS1toS2 = UnitVec3(from1to2/c.distance, true);
        // TODO: check for zero distance and barf
    }
    void realizeMotion       (const State&) const { }
    void realizeDynamics     (const State&) const { }
    void realizeReaction     (const State&) const { }

private:
    const Parameters& getStateVariables(const State& s) const {
        // TODO: assert(parametrized), valid stateIndex
        // TODO: stage must be >= Modeled
        return Value<Parameters>::downcast(s.getDiscreteVariable(stateIndex));
    }

    const Results& getCacheEntries(const State& s) const {
        // TODO: assert valid cacheIndex
        // TODO: stage >= Configured
        return Value<Results>::downcast(s.getCacheEntry(cacheIndex));
    }

    Parameters& updStateVariables(State& s) const {
        // TODO: assert(parametrized), valid stateIndex
        // TODO: stage must be >= Modeled
        return Value<Parameters>::downcast(s.updDiscreteVariable(stateIndex));
    }

    Results& updCacheEntries(const State& s) const {
        // TODO: assert valid cacheIndex
        // TODO: stage >= Time (i.e., config in progress)
        return Value<Results>::downcast(s.updCacheEntry(cacheIndex));
    }

private:
    // If parametrized, we'll store one of these in the State
    class Parameters {
    public:
        Parameters(const Vec3& s1, const Vec3& s2,
                       const Real& len, const Real& k)
          : station1(s1), station2(s2), naturalLength(len), stiffness(k) { }

        String toString() const {
            std::ostringstream s; 
            s << "<Parameters>\n"; 
            s << "<Vec3 name=station1>" << station1 << "</Vec3>\n";
            s << "<Vec3 name=station2>" << station2 << "</Vec3>\n";
            s << "<Real name=naturalLEngth>" << naturalLength << "</Real>\n";
            s << "<Real name=stiffness>" << stiffness << "</Real>\n";
            s << "</Parameters>\n";
            return s.str(); 
        }

        Vec3 station1;
        Vec3 station2;
        Real naturalLength;
        Real stiffness;
    };

    int            body1, body2;
    Parameters defaultState;

    bool      parametrized;
    int       stateIndex;
    int       cacheIndex;

    class Results {
    public:
        Results() : distance(CNT<Real>::getNaN()), fromS1toS2() { }

        String toString() const {
            std::ostringstream s; 
            s << "<Results>\n"; 
            s << "<Real name=distance>" << distance << "</Real>\n";
            s << "<UnitVec3 name=fromS1toS2>" << fromS1toS2 << "</UnitVec3>\n";
            s << "</Results>\n";
            return s.str(); 
        }

        Real     distance;
        UnitVec3 fromS1toS2;
    };

    friend std::ostream& operator<<(std::ostream&, const Parameters& s);
    friend std::ostream& operator<<(std::ostream&, const Results& s);
};

// TODO: these are required to use Value<T>
inline std::ostream& 
operator<<(std::ostream& o, const LinearSpring::Parameters& s) {return o<<s.toString();}
inline std::ostream& 
operator<<(std::ostream& o, const LinearSpring::Results& c) {return o<<c.toString();}

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_H_
