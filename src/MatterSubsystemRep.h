#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_

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
 * Define the private implementation MatterSubsystemRep of a MatterSubsystem, 
 * a still-abstract class derived from abstract base class SubsystemRep.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SubsystemRep.h"

namespace SimTK {

class State;

class MatterSubsystemRep : public SubsystemRep {
public:
    MatterSubsystemRep(const String& name, const String& version)
      : SubsystemRep(name,version), forceSubsys(-1)
    {
    }
    virtual ~MatterSubsystemRep() { }

    void setForceSubsystemIndex(int subsys) {
        assert(subsys >= 0);
        assert(forceSubsys == -1);
        forceSubsys = subsys;
    }
    int getForceSubsystemIndex() const {
        assert(forceSubsys >= 0);
        return forceSubsys;
    }

    const ForceSubsystem& getForceSubsystem() const {
        return ForceSubsystem::downcast(getSystem().getSubsystem(forceSubsys));
    }

    // Topological information.
    virtual int getNBodies()      const = 0;    // includes ground, also # tree joints+1
    virtual int getNMobilities()  const = 0;
    virtual int getNConstraints() const = 0;    // i.e., constraint elements (multiple equations)

    virtual int         getParent  (int bodyNum)           const = 0;
    virtual Array<int>  getChildren(int bodyNum)           const = 0;

    virtual const Transform&  getJointFrame(const State&, int bodyNum) const = 0;
    virtual const Transform&  getJointFrameOnParent(const State&, int bodyNum) const = 0;

    virtual const Vec3&       getBodyCenterOfMass (const State&, int bodyNum) const = 0;
    virtual const Transform&  getBodyConfiguration(const State&, int bodyNum) const = 0;
    virtual const SpatialVec& getBodyVelocity     (const State&, int bodyNum) const = 0;

    // These are simple operators for helping force subsystems put their forces in the 
    // right slots.
    virtual void addInGravity(const State& s, const Vec3& g, Vector_<SpatialVec>& rigidBodyForces) const = 0;
    virtual void addInPointForce(const State& s, int body, const Vec3& stationInB, const Vec3& forceInG,
                                 Vector_<SpatialVec>& rigidBodyForces) const = 0;
    virtual void addInBodyTorque(const State& s, int body, const Vec3& torqueInG, 
                                 Vector_<SpatialVec>& rigidBodyForces) const = 0;
    virtual void addInMobilityForce(const State& s, int body, int axis, const Real& r, 
                                    Vector& mobilityForces) const = 0;  
    
    virtual const Real& getJointQ(const State&, int body, int axis) const = 0;
    virtual const Real& getJointU(const State&, int body, int axis) const = 0;

    virtual void setJointQ(State&, int body, int axis, const Real&) const = 0;
    virtual void setJointU(State&, int body, int axis, const Real&) const = 0;

    virtual const Transform& getMobilizerConfiguration(const State&, int body) const = 0;
    virtual const SpatialVec& getMobilizerVelocity(const State&, int body) const = 0;
    virtual void setMobilizerConfiguration(State&, int body, const Transform& X_JbJ) const = 0;
    virtual void setMobilizerVelocity(State&, int body, const SpatialVec& V_JbJ) const = 0;

    virtual const Vector& getQConstraintErrors(const State&) const {

        static Vector dummy;
        return dummy;
    }
    virtual Real calcQConstraintNorm(const State&) const {
        return 0;
    }
    virtual const Vector& getUConstraintErrors(const State&) const {
        static Vector dummy;
        return dummy;
    }
    virtual Real calcUConstraintNorm(const State&) const {
        return 0;
    }
    virtual const Vector& getUDotConstraintErrors(const State&) const {
        static Vector dummy;
        return dummy;
    }
    virtual Real calcUDotConstraintNorm(const State&) const {
        return 0;
    }
    virtual bool projectQConstraints(State&, Vector& y_err, Real tol, Real targetTol) const {
        return false;
    }
    virtual bool projectUConstraints(State&, Vector& y_err, Real tol, Real targetTol) const {
        return false;
    }

    SimTK_DOWNCAST(MatterSubsystemRep, SubsystemRep);
private:
    int forceSubsys;
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

    SimTK_DOWNCAST(VisualizationSubsystemRep, SubsystemRep);

private:
    // per-body decoration lists
    Array< Array<DecorativeGeometry> > decorations;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_REP_H_
