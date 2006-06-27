#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * Define the private implementation of the MultibodySystem
 * class (a kind of System).
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/MatterSubsystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SystemRep.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

/**
 * This is a global variable for the MultibodySystem, stored as a cache
 * entry in the system's State. Technically it is owned by the
 * MultibodySystemGlobalSubsystem, but it is manipulated directly by
 * the MultibodySystem.
 */
struct DynamicsCache {
    explicit DynamicsCache(int n) // n is # matter subsystems
      : rigidBodyForces(n), particleForces(n), mobilityForces(n),
        potentialEnergy(CNT<Real>::getNaN())
    { }
    std::vector< Vector_<SpatialVec> > rigidBodyForces;
    std::vector< Vector_<Vec3> >       particleForces;
    std::vector< Vector >              mobilityForces;
    Real                               potentialEnergy;
};
// Useless, but required by Value<T>.
inline std::ostream& operator<<(std::ostream& o, const DynamicsCache&) 
{assert(false);return o;}

/**
 * This is the subsystem used by a MultibodySystem to manage global state
 * calculations like forces and potential energy.
 */
class MultibodySystemGlobalSubsystemRep : public SubsystemRep {
    // Topological variables

    mutable int dynamicsCacheIndex;         // where in state to find our stuff
    mutable bool built;

    const DynamicsCache& getDynamicsCache(const State& s) const {
        assert(built);
        return Value<DynamicsCache>::downcast(
            getCacheEntry(s,dynamicsCacheIndex)).get();
    }
    DynamicsCache& updDynamicsCache(const State& s) const {
        assert(built);
        return Value<DynamicsCache>::downcast(
            updCacheEntry(s,dynamicsCacheIndex)).upd();
    }
public:
    MultibodySystemGlobalSubsystemRep()
      : SubsystemRep("MultibodySystemGlobalSubsystem", "0.0.1"),
        dynamicsCacheIndex(-1), built(false)
    {
    }

    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    const Vector_<SpatialVec>& getRigidBodyForces(const State& s, int i) const {
        return getDynamicsCache(s).rigidBodyForces[i];
    }
    const Vector_<Vec3>& getParticleForces(const State& s, int i) const {
        return getDynamicsCache(s).particleForces[i];
    }
    const Vector& getMobilityForces(const State& s, int i) const {
        return getDynamicsCache(s).mobilityForces[i];
    }
    const Real& getPotentialEnergy(const State& s) const {
        return getDynamicsCache(s).potentialEnergy;
    }
    
    Vector_<SpatialVec>& updRigidBodyForces(const State& s, int i) const {
        return updDynamicsCache(s).rigidBodyForces[i];
    }
    Vector_<Vec3>& updParticleForces(const State& s, int i) const {
        return updDynamicsCache(s).particleForces[i];
    }
    Vector& updMobilityForces(const State& s, int i) const {
        return updDynamicsCache(s).mobilityForces[i];
    }
    Real& updPotentialEnergy(const State& s) const {
        return updDynamicsCache(s).potentialEnergy;
    }

    MultibodySystemGlobalSubsystemRep* cloneSubsystemRep() const {
        return new MultibodySystemGlobalSubsystemRep(*this);
    }


    void realizeConstruction(State& s) const {
        const MultibodySystem& mbs = getMultibodySystem();
        dynamicsCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), 
            Stage::Dynamics, 
            new Value<DynamicsCache>(DynamicsCache(mbs.getNMatterSubsystems())));

        built = true;
    }

    // At dynamics stage we make sure the force arrays are all the right length
    // and then initialize them and potential energy to zero. We expect the
    // force subsystems to then fill these in as they are realized to this stage.
    void realizeDynamics(const State& s) const {
        DynamicsCache& dc = updDynamicsCache(s);
        dc.potentialEnergy = 0;
        const MultibodySystem& mbs = getMultibodySystem();
        for (int i=0; i < mbs.getNMatterSubsystems(); ++i) {
            const MatterSubsystem& matter = mbs.getMatterSubsystem(i);
            dc.rigidBodyForces[i].resize(matter.getNBodies());
            dc.particleForces[i].resize(matter.getNParticles());
            dc.mobilityForces[i].resize(matter.getNMobilities());
            dc.rigidBodyForces[i] = SpatialVec(Vec3(0), Vec3(0));
            dc.particleForces[i] = Vec3(0);
            dc.mobilityForces[i] = 0;
        }
    }


    // no need for other realize() methods
    SimTK_DOWNCAST(MultibodySystemGlobalSubsystemRep, SubsystemRep);
};

class MultibodySystemGlobalSubsystem : public Subsystem {
public:
    MultibodySystemGlobalSubsystem() : Subsystem() {
        rep = new MultibodySystemGlobalSubsystemRep();
        rep->setMyHandle(*this);
    }

    SimTK_PIMPL_DOWNCAST(MultibodySystemGlobalSubsystem, Subsystem);
    const MultibodySystemGlobalSubsystemRep& getRep() const;
    MultibodySystemGlobalSubsystemRep&       updRep();
};


/**
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MatterSubsystem and a ForceSubsystem.
 */
class MultibodySystemRep : public SystemRep {
    enum {
        GlobalSubsystemIndex            = 0,
        MatterSubsystemIndex            = 1,
        ForceSubsystemIndex             = 2,
        AnalyticGeometrySubsystemIndex  = 3,
        MassPropertiesSubsystemIndex    = 4,
        VisualizationSubsystemIndex     = 5
    };
public:
    MultibodySystemRep()
      : SystemRep(3, "MultibodySystem", "0.0.1")  // TODO: should be 6
    {
        globalSubs = GlobalSubsystemIndex;
        //TODO: just one each for now but that has to change, esp. for forces
        matterSubs.push_back(MatterSubsystemIndex);
        forceSubs.push_back(ForceSubsystemIndex);
    }
    ~MultibodySystemRep() {
    }

    
    MatterSubsystem& addMatterSubsystem(MatterSubsystem& m) {
        // TODO: allow more than one
        Subsystem& s = takeOverSubsystem(MatterSubsystemIndex, m);
        return MatterSubsystem::updDowncast(s);
    }
    ForceSubsystem& addForceSubsystem(ForceSubsystem& f) {
        // TODO: allow more than one
        Subsystem& s = takeOverSubsystem(ForceSubsystemIndex, f);
        return ForceSubsystem::updDowncast(s);
    }

    int getNMatterSubsystems() const {return matterSubs.size();}
    int getNForceSubsystems()  const {return forceSubs.size();}
    const MatterSubsystem& getMatterSubsystem(int i) const {
        return MatterSubsystem::downcast(getSubsystem(matterSubs[i]));
    }
    const ForceSubsystem& getForceSubsystem(int i) const {
        return ForceSubsystem::downcast(getSubsystem(forceSubs[i]));
    }
    const MultibodySystemGlobalSubsystem& getGlobalSubsystem() const {
        return MultibodySystemGlobalSubsystem::downcast(getSubsystem(globalSubs));
    }

    MatterSubsystem& updMatterSubsystem(int i) {
        return MatterSubsystem::updDowncast(updSubsystem(matterSubs[i]));
    }
    ForceSubsystem& updForceSubsystem(int i) {
        return ForceSubsystem::updDowncast(updSubsystem(forceSubs[i]));
    }
    MultibodySystemGlobalSubsystem& updGlobalSubsystem() {
        return MultibodySystemGlobalSubsystem::updDowncast(updSubsystem(globalSubs));
    }
    // Global state variables dealing with interaction between forces & matter

    // Responses available when the global subsystem is advanced to Dynamics stage.
    const Vector_<SpatialVec>& getRigidBodyForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().getRigidBodyForces(s,i);
    }
    const Vector_<Vec3>& getParticleForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().getParticleForces(s,i);
    }
    const Vector& getMobilityForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().getMobilityForces(s,i);
    }
    const Real& getPotentialEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().getPotentialEnergy(s);
    }

    // Dynamics stage cache entries of the global subsystem. Accessing these drops
    // the global stage to Dynamics-1, i.e. Moving.
    Vector_<SpatialVec>& updRigidBodyForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().updRigidBodyForces(s,i);
    }
    Vector_<Vec3>& updParticleForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().updParticleForces(s,i);
    }
    Vector& updMobilityForces(const State& s, int i) const {
        return getGlobalSubsystem().getRep().updMobilityForces(s,i);
    }
    Real& updPotentialEnergy(const State& s) const {
        return getGlobalSubsystem().getRep().updPotentialEnergy(s);
    }

    // Each matter subsystem's constraints are independent of the other
    // matter subsystems, so we can project them individually. Do all
    // "q" projections before any "u" projections, though.
    bool project(State& s, Vector& y_err,
                 const Real& tol, const Real& dontProjectFac, 
                 const Real& targetTol) const 
    {
        bool anyChange = false;

        realize(s, Stage::Timed);

        for (int i=0; i < getNMatterSubsystems(); ++i) {
            const MatterSubsystem& mech = getMatterSubsystem(i);
            mech.realize(s, Stage::Configured);
            const Real qerr = mech.calcQConstraintNorm(s);
            if (dontProjectFac==0 || qerr > tol*dontProjectFac) {
                if (mech.projectQConstraints(s, y_err, tol, targetTol))
                    anyChange = true;
            }
        }

        realize(s, Stage::Configured);

        for (int i=0; i < getNMatterSubsystems(); ++i) {
            const MatterSubsystem& mech = getMatterSubsystem(i);
            mech.realize(s, Stage::Moving);
            const Real uerr = mech.calcUConstraintNorm(s);
            if (dontProjectFac==0 || uerr > tol*dontProjectFac) {
                if (mech.projectUConstraints(s, y_err, tol, targetTol))
                    anyChange = true;
            }
        }

        return anyChange;
    }



    // pure virtual
    MultibodySystemRep* cloneSystemRep() const {return new MultibodySystemRep(*this);}


    void realizeConstruction(State& s) const {
        MultibodySystemRep& mutableThis = *const_cast<MultibodySystemRep*>(this);
        // Create the global subsystem
        (void)mutableThis.takeOverSubsystem(GlobalSubsystemIndex, MultibodySystemGlobalSubsystem());

        // Help the subsystems find each other. TODO delete
        mutableThis.updMatterSubsystem(0).setForceSubsystemIndex( ForceSubsystemIndex );
        mutableThis.updForceSubsystem(0).setMatterSubsystemIndex( MatterSubsystemIndex );

        getGlobalSubsystem().realize(s, Stage::Built);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Built);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Built);
    }
    void realizeModeling(State& s) const {
        getGlobalSubsystem().realize(s, Stage::Modeled);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Modeled);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Modeled);
    }
    void realizeParameters(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Parametrized);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Parametrized);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Parametrized);
    }
    void realizeTime(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Timed);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Timed);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Timed);
    }
    void realizeConfiguration(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Configured);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Configured);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Configured);
    }
    void realizeMotion(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Moving);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Moving);
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Moving);
    }
    void realizeDynamics(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Dynamics);
        // note order: forces first (TODO: does that matter?)
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Dynamics);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Dynamics);
    }
    void realizeReaction(const State& s) const {
        getGlobalSubsystem().realize(s, Stage::Reacting);
        // note order: forces first (TODO: does that matter?)
        for (int i=0; i < getNForceSubsystems(); ++i)
            getForceSubsystem(i).realize(s, Stage::Reacting);
        for (int i=0; i < getNMatterSubsystems(); ++i)
            getMatterSubsystem(i).realize(s, Stage::Reacting);
    }



    SimTK_DOWNCAST(MultibodySystemRep, SystemRep);
private:
    int globalSubs;                 // index of global subsystem
    std::vector<int> matterSubs;    // indices of matter subsystems
    std::vector<int> forceSubs;     // indices of force subsystems
};



} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
