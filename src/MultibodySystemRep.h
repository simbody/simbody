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
 * The job of the MultibodySystem class is to coordinate the activities of a
 * MatterSubsystem and a ForceSubsystem.
 */
class MultibodySystemRep : public SystemRep {
    enum {
        SystemSubsystemIndex            = 0,
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
    }
    ~MultibodySystemRep() {
    }

    
    bool project(State& s, Vector& y_err,
                 const Real& tol, const Real& dontProjectFac, 
                 const Real& targetTol) const 
    {
        const MatterSubsystem& mech = getMatterSubsystem();
        bool anyChange = false;

        realize(s, Stage::Timed);
        mech.realize(s, Stage::Configured);
        const Real qerr = mech.getQConstraintNorm(s);
        if (qerr > tol*dontProjectFac) {
            if (mech.projectQConstraints(s, y_err, tol, targetTol))
                anyChange = true;
        }
        realize(s, Stage::Configured);
        mech.realize(s, Stage::Moving);
        const Real uerr = mech.getUConstraintNorm(s);
        if (uerr > tol*dontProjectFac) {
            if (mech.projectUConstraints(s, y_err, tol, targetTol))
                anyChange = true;
        }

        return anyChange;
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
        // Help the subsystems find each other.
        MultibodySystemRep& mutableThis = *const_cast<MultibodySystemRep*>(this);
        mutableThis.updMatterSubsystem().setForceSubsystemIndex( ForceSubsystemIndex );
        mutableThis.updForceSubsystem().setMatterSubsystemIndex( MatterSubsystemIndex );

        getMatterSubsystem().realize(s, Stage::Built);
        getForceSubsystem().realize(s, Stage::Built);
    }
    void realizeModeling(State& s) const {
        getMatterSubsystem().realize(s, Stage::Modeled);
        getForceSubsystem().realize(s, Stage::Modeled);
    }
    void realizeParameters(const State& s) const {
        getMatterSubsystem().realize(s, Stage::Parametrized);
        getForceSubsystem().realize(s, Stage::Parametrized);
    }
    void realizeTime(const State& s) const {
        getMatterSubsystem().realize(s, Stage::Timed);
        getForceSubsystem().realize(s, Stage::Timed);
    }
    void realizeConfiguration(const State& s) const {
        getMatterSubsystem().realize(s, Stage::Configured);
        getForceSubsystem().realize(s, Stage::Configured);
    }
    void realizeMotion(const State& s) const {
        getMatterSubsystem().realize(s, Stage::Moving);
        getForceSubsystem().realize(s, Stage::Moving);
    }
    void realizeDynamics(const State& s) const {
        getForceSubsystem().realize(s, Stage::Dynamics); // note order
        getMatterSubsystem().realize(s, Stage::Dynamics);
    }
    void realizeReaction(const State& s) const {
        getForceSubsystem().realize(s, Stage::Reacting);
        getMatterSubsystem().realize(s, Stage::Reacting);
    }

    MatterSubsystem& setMatterSubsystem(MatterSubsystem& m) {
        bodies.resize(m.getNBodies());
        Subsystem& s = takeOverSubsystem(MatterSubsystemIndex, m);
        return MatterSubsystem::updDowncast(s);
    }
    ForceSubsystem& setForceSubsystem(ForceSubsystem& f) {
        Subsystem& s = takeOverSubsystem(ForceSubsystemIndex, f);
        return ForceSubsystem::updDowncast(s);
    }

    const MatterSubsystem& getMatterSubsystem() const {
        return MatterSubsystem::downcast(getSubsystem(MatterSubsystemIndex));
    }
    const ForceSubsystem& getForceSubsystem() const {
        return ForceSubsystem::downcast(getSubsystem(ForceSubsystemIndex));
    }


    MatterSubsystem& updMatterSubsystem() {
        return MatterSubsystem::updDowncast(updSubsystem(MatterSubsystemIndex));
    }
    ForceSubsystem& updForceSubsystem() {
        return ForceSubsystem::updDowncast(updSubsystem(ForceSubsystemIndex));
    }

    const Array<AnalyticGeometry>&   getBodyAnalyticGeometry(int bodyNum) const {
        return bodies[bodyNum].aGeom;
    }
    const Array<DecorativeGeometry>& getBodyDecorativeGeometry(int bodyNum) const {
        return bodies[bodyNum].dGeom;
    }

    SimTK_DOWNCAST(MultibodySystemRep, SystemRep);
private:
    struct PerBodyInfo {
        Array<AnalyticGeometry>   aGeom;
        Array<DecorativeGeometry> dGeom;
    };
    Array<PerBodyInfo> bodies;
};



} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_REP_H_
