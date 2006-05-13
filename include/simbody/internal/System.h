#ifndef SimTK_SIMBODY_SYSTEM_H_
#define SimTK_SIMBODY_SYSTEM_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

// Stubs.

class SimTK_SIMBODY_API System {
public:
    System() { }
    virtual ~System() { }

    virtual System* cloneSystem() const = 0;

    virtual void realizeConstruction (State& s)       const = 0;
    virtual void realizeModeling     (State& s)       const = 0;
    virtual void realizeParameters   (const State& s) const { }
    virtual void realizeTime         (const State& s) const { }
    virtual void realizeConfiguration(const State& s) const { }
    virtual void realizeMotion       (const State& s) const { }
    virtual void realizeDynamics     (const State& s) const { }
    virtual void realizeReaction     (const State& s) const { }

    void realize(const State& s, Stage g) const {
        while (s.getStage() < g) {
            switch (s.getStage()) {
            case Stage::Allocated:    realizeConstruction(const_cast<State&>(s)); break;
            case Stage::Built:        realizeModeling    (const_cast<State&>(s)); break;
            case Stage::Modeled:      realizeParameters(s);    break;
            case Stage::Parametrized: realizeTime(s);          break;
            case Stage::Timed:        realizeConfiguration(s); break;
            case Stage::Configured:   realizeMotion(s);        break;
            case Stage::Moving:       realizeDynamics(s);      break;
            case Stage::Dynamics:     realizeReaction(s);      break;
            default: assert(!"System::realize(): bad stage");
            }
            s.advanceToStage(s.getStage().next());
        }
    }

private:

};

class SimTK_SIMBODY_API Study {
public:
    Study(const System& sys)
      : system(sys.cloneSystem())
    {
        system->realizeConstruction(state);
    }

    ~Study() {
        delete system;
    }

    const System& getSystem() const {return *system;}
    const State&  getState()  const {return state;}
    State&        updState()        {return state;}
private:
    System* system;
    State   state;
};



class SimTK_SIMBODY_API Subsystem {
public:
    virtual ~Subsystem() { }

    virtual Subsystem* cloneSubsystem() const = 0;
    virtual void endConstruction() { }
    virtual void realizeConstruction(State&) const = 0;
    virtual void realizeModeling(State&) const = 0;

private:
};

class MechanicalForcesSubsystem;
class SimTK_SIMBODY_API MechanicalSubsystem : public Subsystem {
public:
    MechanicalSubsystem() { }
    virtual ~MechanicalSubsystem() { }

    Subsystem* cloneSubsystem() const {return cloneMechanicalSubsystem();}

    virtual MechanicalSubsystem* cloneMechanicalSubsystem() const = 0;

    virtual int getNBodies() const = 0;
    virtual const Transform& getBodyConfiguration(const State&, int) const = 0;


    virtual void realizeParameters   (const State&) const { }
    virtual void realizeTime         (const State&) const { }
    virtual void realizeConfiguration(const State&) const { }
    virtual void realizeMotion       (const State&) const { }
    virtual void realizeDynamics     (const State&, const MechanicalForcesSubsystem&) const { }
    virtual void realizeReaction     (const State&, const MechanicalForcesSubsystem&) const { }

    virtual const Real& getJointQ(const State&, int body, int axis) const = 0;
    virtual const Real& getJointU(const State&, int body, int axis) const = 0;

    virtual void setJointQ(State&, int body, int axis, const Real&) const = 0;
    virtual void setJointU(State&, int body, int axis, const Real&) const = 0;

    SimTK_DOWNCAST(MechanicalSubsystem, Subsystem);
private:
};

class SimTK_SIMBODY_API MechanicalForcesSubsystem : public Subsystem {
public:
    MechanicalForcesSubsystem(const MechanicalSubsystem& m) 
        : mech(m) 
    {
    }

    Subsystem* cloneSubsystem() const {return cloneMechanicalForcesSubsystem();}

    virtual MechanicalForcesSubsystem* cloneMechanicalForcesSubsystem() const = 0;

    virtual void realizeParameters   (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeTime         (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeConfiguration(const State&, const MechanicalSubsystem&) const { }
    virtual void realizeMotion       (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeDynamics     (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeReaction     (const State&, const MechanicalSubsystem&) const { }

    SimTK_DOWNCAST(MechanicalForcesSubsystem, Subsystem);
private:
    const MechanicalSubsystem& mech;
};

/**
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MechanicalSubsystem and a MechanicalForcesSubsystem.
 */
class SimTK_SIMBODY_API MultibodySystem : public System {
public:
    MultibodySystem(const MechanicalSubsystem& m, const MechanicalForcesSubsystem& f)
        : mech(m.cloneMechanicalSubsystem()), forces(f.cloneMechanicalForcesSubsystem())
    {
        mech->endConstruction();     // in case user forgot ...
        forces->endConstruction();
    }
    ~MultibodySystem() {
        delete mech; mech=0;
        delete forces; forces=0;
    }

    // TODO: camera facing, screen fixed, calculated geometry (e.g. line between stations
    // on two different bodies, marker at system COM)
    void addAnalyticGeometry  (int bodyNum, const Transform&, const AnalyticGeometry&);
    void addDecorativeGeometry(int bodyNum, const Transform&, Real scale, const DecorativeGeometry&);


    System* cloneSystem() const {return new MultibodySystem(*this);}

    void realizeConstruction(State& s) const {
        mech->realizeConstruction(s);
        forces->realizeConstruction(s);
    }
    void realizeModeling(State& s) const {
        mech->realizeModeling(s);
        forces->realizeModeling(s);
    }
    void realizeParameters(const State& s) const {
        mech->realizeParameters(s);
        forces->realizeParameters(s, *mech);
    }
    void realizeTime(const State& s) const {
        mech->realizeTime(s);
        forces->realizeTime(s, *mech);
    }
    void realizeConfiguration(const State& s) const {
        mech->realizeConfiguration(s);
        forces->realizeConfiguration(s, *mech);
    }
    void realizeMotion(const State& s) const {
        mech->realizeMotion(s);
        forces->realizeMotion(s, *mech);
    }
    void realizeDynamics(const State& s) const {
        forces->realizeDynamics(s, *mech); // note order
        mech->realizeDynamics(s, *forces);
    }
    void realizeReaction(const State& s) const {
        forces->realizeReaction(s, *mech);
        mech->realizeReaction(s, *forces);
    }

    const MechanicalSubsystem&       getMechanicalSubsystem()       const {return *mech;}
    const MechanicalForcesSubsystem& getMechanicalForcesSubsystem() const {return *forces;}

    SimTK_DOWNCAST(MultibodySystem, System);
private:
    MechanicalSubsystem*       mech;
    MechanicalForcesSubsystem* forces;
};

class SimTK_SIMBODY_API MultibodyDynamicsStudy : public Study {
public:
    MultibodyDynamicsStudy(const MultibodySystem& sys)
        : Study(sys)
    {
    }

    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    void advanceTimeBy(const Real& h) { 
        printf("advanceTimeBy(%g) ... TODO!\n", h);
    }
private:
};


} // namespace SimTK

#endif // SimTK_SIMBODY_SYSTEM_H_
