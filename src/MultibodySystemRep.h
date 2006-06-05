#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_

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
 * Define the private implementation of the MultibodySystem
 * class (a kind of System).
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SystemRep.h"

#include <vector>

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;


class MechanicalForcesSubsystemRep;
class MechanicalSubsystemRep : public SubsystemRep {
public:
    MechanicalSubsystemRep() { }
    virtual ~MechanicalSubsystemRep() { }

    // Topological information.
    virtual int getNBodies()      const = 0;    // includes ground, also # tree joints+1
    virtual int getNConstraints() const = 0;    // i.e., constraint elements (multiple equations)

    virtual int         getParent  (int bodyNum)           const = 0;
    virtual Array<int>  getChildren(int bodyNum)           const = 0;

    virtual const Transform&  getJointFrame(const State&, int bodyNum) const = 0;
    virtual const Transform&  getJointFrameOnParent(const State&, int bodyNum) const = 0;

    virtual const Vec3&       getBodyCenterOfMass (const State&, int bodyNum) const = 0;
    virtual const Transform&  getBodyConfiguration(const State&, int bodyNum) const = 0;
    virtual const SpatialVec& getBodyVelocity     (const State&, int bodyNum) const = 0;

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

    SimTK_DOWNCAST(MechanicalSubsystemRep, SubsystemRep);
};

class MechanicalForcesSubsystemRep : public SubsystemRep {
public:
    MechanicalForcesSubsystemRep(const String& name, const String& version,
                                 const MechanicalSubsystem& m) 
        : SubsystemRep(name,version), mech(m) 
    {
    }

    const MechanicalSubsystem& getMechanicalSubsystem() const {return mech;}

    virtual void realizeParameters   (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeTime         (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeConfiguration(const State&, const MechanicalSubsystem&) const { }
    virtual void realizeMotion       (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeDynamics     (const State&, const MechanicalSubsystem&) const { }
    virtual void realizeReaction     (const State&, const MechanicalSubsystem&) const { }

    SimTK_DOWNCAST(MechanicalForcesSubsystemRep, SubsystemRep);
private:
    const MechanicalSubsystem& mech;
};



// Concrete for now.
class Body {
public:
    int getBodyNumber() const;
};

class VisualizationSubsystemRep : public SubsystemRep {
public:
    VisualizationSubsystemRep(const String& name, const String& version) 
      : SubsystemRep(name, version)
    {
    }

    void addDecorativeGeometry(const Body& b, const Transform& X_BG, const DecorativeGeometry& g)
    {
        const int bnum = b.getBodyNumber();
        if (decorations.size() <= bnum)
            decorations.resize(bnum+1);
        decorations[bnum].push_back(g);
        decorations[bnum].back().setPlacement(X_BG);
    }

    const Array<DecorativeGeometry>& getBodyDecorativeGeometry(const Body& b) const {
        static const Array<DecorativeGeometry> empty;
        const int bnum = b.getBodyNumber();
        return bnum < decorations.size() ? decorations[bnum] : empty;
    }

    void realizeConstruction (State&)       const { }
    void realizeModeling     (State&)       const { }
    void realizeParameters   (const State&) const { }
    void realizeTime         (const State&) const { }
    void realizeConfiguration(const State&) const { }
    void realizeMotion       (const State&) const { }
    void realizeDynamics     (const State&) const { }
    void realizeReaction     (const State&) const { }

    SimTK_DOWNCAST(VisualizationSubsystemRep, SubsystemRep);

private:
    // per-body decoration lists
    Array< Array<DecorativeGeometry> > decorations;
};

/**
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MechanicalSubsystem and a MechanicalForcesSubsystem.
 */
class MultibodySystemRep : public SystemRep {
    enum {
        SystemSubsystemIndex            = 0,
        MechanicalSubsystemIndex        = 1,
        MechanicalForcesSubsystemIndex  = 2,
        AnalyticGeometrySubsystemIndex  = 3,
        MassPropertiesSubsystemIndex    = 4,
        VisualizationSubsystemIndex     = 5
    };
public:
    MultibodySystemRep(const MechanicalSubsystem& m, const MechanicalForcesSubsystem& f)
        : SystemRep(6, "MultibodySystem", "0.0.1"), mech(m), forces(f)
    {
        //mech.endConstruction();     // in case user forgot ...
        //forces.endConstruction();

        bodies.resize(mech.getNBodies());
    }
    ~MultibodySystemRep() {
    }

    // pure virtual
    SystemRep* cloneSystemRep() const {return new MultibodySystemRep(*this);}


    // TODO: camera facing, screen fixed, calculated geometry (e.g. line between stations
    // on two different bodies, marker at system COM)
    void addAnalyticGeometry  (int bodyNum, const Transform& X_BG, const AnalyticGeometry& g)
    {
        assert(0 <= bodyNum && bodyNum < bodies.size());
        bodies[bodyNum].aGeom.push_back(g);
        bodies[bodyNum].aGeom.back().setPlacement(X_BG);
    }
    void addDecorativeGeometry(int bodyNum, const Transform& X_BG, const DecorativeGeometry& g)
    {
        assert(0 <= bodyNum && bodyNum < bodies.size());
        bodies[bodyNum].dGeom.push_back(g);
        bodies[bodyNum].dGeom.back().setPlacement(X_BG);
    }

    void realizeConstruction(State& s) const {
        mech.realize(s, Stage::Built);
        forces.realize(s, Stage::Built);
    }
    void realizeModeling(State& s) const {
        mech.realize(s, Stage::Modeled);
        forces.realize(s, Stage::Modeled);
    }
    void realizeParameters(const State& s) const {
        mech.realize(s, Stage::Parametrized);
        forces.realize(s, Stage::Parametrized);
    }
    void realizeTime(const State& s) const {
        mech.realize(s, Stage::Timed);
        forces.realize(s, Stage::Timed);
    }
    void realizeConfiguration(const State& s) const {
        mech.realize(s, Stage::Configured);
        forces.realize(s, Stage::Configured);
    }
    void realizeMotion(const State& s) const {
        mech.realize(s, Stage::Moving);
        forces.realize(s, Stage::Moving);
    }
    void realizeDynamics(const State& s) const {
        forces.realize(s, Stage::Dynamics); // note order
        mech.realize(s, Stage::Dynamics);
    }
    void realizeReaction(const State& s) const {
        forces.realize(s, Stage::Reacting);
        mech.realize(s, Stage::Reacting);
    }

    const MechanicalSubsystem&       getMechanicalSubsystem()       const {return mech;}
    const MechanicalForcesSubsystem& getMechanicalForcesSubsystem() const {return forces;}
    const Array<AnalyticGeometry>&   getBodyAnalyticGeometry(int bodyNum) const {
        return bodies[bodyNum].aGeom;
    }
    const Array<DecorativeGeometry>& getBodyDecorativeGeometry(int bodyNum) const {
        return bodies[bodyNum].dGeom;
    }

    SimTK_DOWNCAST(MultibodySystemRep, SystemRep);
private:
    const MechanicalSubsystem&       mech;
    const MechanicalForcesSubsystem& forces;

    struct PerBodyInfo {
        Array<AnalyticGeometry>   aGeom;
        Array<DecorativeGeometry> dGeom;
    };
    Array<PerBodyInfo> bodies;
};



} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
