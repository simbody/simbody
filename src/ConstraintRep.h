#ifndef SimTK_SIMBODY_CONSTRAINT_REP_H_
#define SimTK_SIMBODY_CONSTRAINT_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**@file
 *
 * Private implementation of Constraint and its included subclasses which
 * represent the built-in constraint types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

class SimTK::SimbodyMatterSubsystem;

class ConstraintNode;
class SimbodyMatterSubsystemRep;

namespace SimTK {

    /////////////////////
    // CONSTRAINT REPS //
    /////////////////////

class Constraint::ConstraintRep {
public:
    ConstraintRep() : myHandle(0), myMatterSubsystemRep(0), myConstraintNode(0) {
    }
    virtual ~ConstraintRep();
    virtual ConstraintRep* clone() const = 0;

    // This creates a constraint node using the appropriate constraint type.
    virtual ConstraintNode* createConstraintNode() const = 0; 

    void realizeTopology(int& nxtQErr, int& nxtUErr, int& nxtMult) const;

    virtual void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
    {
    }

    void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
    {
        // Let the individual constraint deal with any complicated stuff.
        calcDecorativeGeometryAndAppendImpl(s,stage,geom);
    }

    void invalidateTopologyCache() const;
    const ConstraintNode& getMyConstraintNode() const;

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              ConstraintId id);

    bool isInSubsystem() const {
        return myMatterSubsystemRep != 0;
    }

    const SimbodyMatterSubsystemRep& getMyMatterSubsystemRep() const {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }
    SimbodyMatterSubsystemRep& updMyMatterSubsystemRep() {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }

    void setMyHandle(Constraint& h) {myHandle = &h;}
    const Constraint& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class Constraint;
    Constraint* myHandle;	// the owner handle of this rep

        // TOPOLOGY "STATE"

    // These data members are filled in once the Constraint is added to
    // a MatterSubsystem.
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    ConstraintId               myConstraintId; // id within the subsystem

        // TOPOLOGY "CACHE"
    // A ConstraintNode is created for each Constraint during realizeTopology().
    // Think of it as the *computational* form of the Constraint; the Constraint
    // itself exists to make a nice API for the programmer. In fact is is possible
    // for several different Constraint classes to use the same underling computational
    // form while providing a different API.
    //
    // This is a pointer to an abstract object whose heap space is *owned* by
    // this Constraint. Be sure to delete it upon destruction and whenever
    // topology is re-realized.
    mutable ConstraintNode* myConstraintNode;
};


class Constraint::Rod::RodRep : public Constraint::ConstraintRep {
public:
    RodRep() : defaultPoint1(0), defaultPoint2(0), defaultRodLength(1) { }
    RodRep* clone() const { return new RodRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(RodRep, ConstraintRep);
private:
    friend class Constraint::Rod;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRodLength;
};


class Constraint::PointInPlane::PointInPlaneRep : public Constraint::ConstraintRep {
public:
    PointInPlaneRep()
      : defaultPlaneNormal(), defaultPlaneHeight(0), defaultFollowerPoint(0),
        planeHalfWidth(1), pointRadius(0.05) 
    { }
    PointInPlaneRep* clone() const { return new PointInPlaneRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const;

    void setPlaneDisplayHalfWidth(Real h) {
        // h <= 0 means don't display plane
        invalidateTopologyCache();
        planeHalfWidth = h > 0 ? h : 0;
    }
    Real getPlaneDisplayHalfWidth() const {return planeHalfWidth;}

    void setPointDisplayRadius(Real r) {
        // r <= 0 means don't display point
        invalidateTopologyCache();
        pointRadius= r > 0 ? r : 0;
    }
    Real getPointDisplayRadius() const {return pointRadius;}

    SimTK_DOWNCAST(PointInPlaneRep, ConstraintRep);
private:
    friend class Constraint::PointInPlane;

    MobilizedBodyId planeBody;    // B1
    MobilizedBodyId followerBody; // B2
    UnitVec3        defaultPlaneNormal; // on body 1, exp. in B1 frame
    Real            defaultPlaneHeight;
    Vec3            defaultFollowerPoint; // on body 2, exp. in B2 frame

    // These are just for visualization
    Real planeHalfWidth;
    Real pointRadius;
};


class Constraint::Ball::BallRep : public Constraint::ConstraintRep {
public:
    BallRep() : defaultPoint1(0), defaultPoint2(0), defaultRadius(0.1) { }
    BallRep* clone() const { return new BallRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const;

    void setDefaultRadius(Real r) {
        // r <= 0 means don't display
        invalidateTopologyCache();
        defaultRadius = r > 0 ? r : 0;
    }
    Real getDefaultRadius() const {return defaultRadius;}

    SimTK_DOWNCAST(BallRep, ConstraintRep);
private:
    friend class Constraint::Ball;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRadius; // used for visualization only
};



class Constraint::Weld::WeldRep : public Constraint::ConstraintRep {
public:
    WeldRep() { } // transforms are Identity
    WeldRep* clone() const { return new WeldRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(WeldRep, ConstraintRep);
private:
    friend class Constraint::Weld;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Transform       defaultFrame1; // on body 1, relative to B1 frame
    Transform       defaultFrame2; // on body 2, relative to B2 frame};
};


class Constraint::Custom::CustomRep : public Constraint::ConstraintRep {
public:
    CustomRep* clone() const { return new CustomRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(CustomRep, ConstraintRep);
private:
    friend class Constraint::Custom;

    // TODO: God only knows what goes here!
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_REP_H_



