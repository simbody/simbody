#ifndef SimTK_MATTER_SUBSYSTEM_REP_H_
#define SimTK_MATTER_SUBSYSTEM_REP_H_

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
 * Define the private implementation MatterSubsystemRep of a MatterSubsystem, 
 * a still-abstract class derived from abstract base class SubsystemRep.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "SubsystemRep.h"

namespace SimTK {

class State;

class MatterSubsystemRep : public SubsystemRep {
public:
    MatterSubsystemRep(const String& name, const String& version)
      : SubsystemRep(name,version)
    {
    }
    virtual ~MatterSubsystemRep() { }

    // Return the MultibodySystem which owns this MatterSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

        // TOPOLOGY STAGE //
    virtual int getNBodies()      const = 0;    // includes ground, also # mobilizers (tree joints) +1
    virtual int getNParticles()   const {return 0;} // TODO
    virtual int getNMobilities()  const = 0;
    virtual int getNConstraints() const = 0;    // i.e., constraint elements (multiple equations)

    virtual BodyId         getParent  (BodyId bodyNum)           const = 0;
    virtual Array<BodyId>  getChildren(BodyId bodyNum)           const = 0;

        // MODEL STAGE //

    // These report the start index and number of generalized coordinates q or generalized
    // speeds u associated with a body's mobilizer. These are indices into this subsystem's
    // Q or U allocation in the state, and can also be used other arrays which have the
    // same dimensions. For example, mobilizer force arrays have the same dimension as
    // the Us (that is, their length is the total mobility of the system).
    virtual void findMobilizerQs(const State& s, BodyId body, int& qStart, int& nq) const = 0;
    virtual void findMobilizerUs(const State& s, BodyId body, int& uStart, int& nu) const = 0;

    virtual void setMobilizerTransform  (State&, BodyId, const Transform& X_MbM) const = 0;
    virtual void setMobilizerRotation   (State&, BodyId, const Rotation&  R_MbM) const = 0;
    virtual void setMobilizerTranslation(State&, BodyId, const Vec3&      T_MbM, 
                                         bool dontChangeOrientation)             const = 0;

    virtual void setMobilizerVelocity       (State&, BodyId, const SpatialVec& V_MbM) const = 0;
    virtual void setMobilizerAngularVelocity(State&, BodyId, const Vec3&       w_MbM) const = 0;
    virtual void setMobilizerLinearVelocity (State&, BodyId, const Vec3&       v_MbM,
                                             bool dontChangeAngularVelocity)          const = 0;


        // INSTANCE STAGE //

    virtual const MassProperties& getBodyMassProperties(const State& s, BodyId body) const = 0;
    virtual const Transform&      getMobilizerFrame(const State&, BodyId) const = 0;
    virtual const Transform&      getMobilizerFrameOnParent(const State&, BodyId) const = 0;

    virtual const Vector&     getParticleMasses(const State&) const { // TODO
        static Vector v;
        return v;
    }

        // POSITION, VELOCITY, ACCELERATION STAGES //
    virtual const Transform&  getBodyTransform(const State&, BodyId) const = 0;
    virtual const SpatialVec& getBodyVelocity(const State&, BodyId) const = 0;
    virtual const SpatialVec& getBodyAcceleration(const State&, BodyId) const = 0; 

    virtual const Transform& getMobilizerTransform(const State&, BodyId) const = 0;
    virtual const SpatialVec& getMobilizerVelocity(const State&, BodyId) const = 0;

    virtual const Vector_<Vec3>& getParticleLocations(const State&) const { // TODO
        static Vector_<Vec3> v;
        return v;
    }

    virtual Real calcQConstraintNorm(const State&) const {
        return 0;
    }
    virtual Real calcUConstraintNorm(const State&) const {
        return 0;
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

#endif // SimTK_MATTER_SUBSYSTEM_REP_H_
