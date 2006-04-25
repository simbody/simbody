#ifndef SimTK_SIMBODY_FORCES_H_
#define SimTK_SIMBODY_FORCES_H_

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"

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

    void realizeConstruction(State&) const { }
    void realizeModeling(State&) const { }

    void setGravity(const Vec3& g) { defaultGravity=g; }

    SimTK_DOWNCAST2(BasicMechanicalForceElements, MechanicalForcesSubsystem, Subsystem);
private:
    Vec3 defaultGravity;
};


/*
 * This class holds a full set of applied forces, including body
 * forces (by which we mean forces and torques), and joint forces..
 *
 * Body forces are expressed in the ground frame. Joints constitute
 * their own frames, and joint forces are applied as per-joint-DOF
 * scalars whose meanings depend on the quirks of the individual
 * joints.
 *
 * Note that this is a low-level container class; it does not know
 * the joint types or current configuration of the bodies. It just
 * knows how many items of each type there are. It must be given to a
 * SimbodyTree, along with the current state, in order to be
 * interpreted correctly.
 *
class SimbodyForceSet {
public:
    SimbodyForceSet() { }
    SimbodyForceSet(int nBody, int nDOFs) : bodyForces(nBody), jointForces(nDOFs) { }
    // default copy, assignment, destructor

    void resize(int nBody, int nDOFs) 
       {bodyForces.resize(nBody); jointForces.resize(nDOFs);}

    /// Set all forces to zero.
    void clear() {clearBodyForces(); clearJointForces();}

    /// Clear only the body forces, leaving joint forces alone.
    void clearBodyForces()  {bodyForces.setToZero();}

    /// Clear only the joint forces, leaving the body forces alone.
    void clearJointForces() {jointForces.setToZero();}

    const Vector_<SpatialVec>& getBodyForces() const {return bodyForces;}
    Vector_<SpatialVec>&       updBodyForces()       {return bodyForces;}

    const Vector& getJointForces() const {return jointForces;}
    Vector&       updJointForces()       {return jointForces;}

    // Named assignment operators; prefer the actual operators in C++
    SimbodyForceSet& assign  (const SimbodyForceSet& f) {return (*this = f);}
    SimbodyForceSet& add     (const SimbodyForceSet& f) 
       {bodyForces+=f.bodyForces; jointForces+=f.jointForces; return *this;}
    SimbodyForceSet& subtract(const SimbodyForceSet& f) 
       {bodyForces-=f.bodyForces; jointForces-=f.jointForces; return *this;}
    SimbodyForceSet& scale   (const Real& s) {bodyForces*=s; jointForces*=s; return *this;}

    SimbodyForceSet& operator+=(const SimbodyForceSet& f) {return add(f);}
    SimbodyForceSet& operator-=(const SimbodyForceSet& f) {return subtract(f);}
    SimbodyForceSet& operator*=(const Real& s)            {return scale(s);}
private:
    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3>       particleForces;
    Vector_<Real>       jointForces;
};
*/
/*
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
*/

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_H_
