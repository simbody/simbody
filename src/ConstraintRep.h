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
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem_Subtree.h"

#include <map>
#include <vector>
#include <utility>  // std::pair

class ConstraintNode;
class SimbodyMatterSubsystemRep;

namespace SimTK {

    /////////////////////
    // CONSTRAINT REPS //
    /////////////////////

class Constraint::ConstraintRep {
public:
    ConstraintRep() : myHandle(0), myMatterSubsystemRep(0), 
        defaultMp(0), defaultMv(0), defaultMa(0),
        myConstraintNode(0) 
    {
    }

    ConstraintRep(int mp, int mv, int ma) : myHandle(0), myMatterSubsystemRep(0), 
        defaultMp(mp), defaultMv(mv), defaultMa(ma),
        myConstraintNode(0) 
    {
    }

    void setDefaultNumConstraints(int mp, int mv, int ma) {
        assert(mp >= 0 && mv >= 0 && ma >= 0);
        invalidateTopologyCache();
        defaultMp = mp;
        defaultMv = mv;
        defaultMa = ma;
    }

    typedef std::map<MobilizedBodyId,ConstrainedBodyId> Mobilized2ConstrainedMap;

    ConstrainedBodyId addConstrainedBody(const MobilizedBody& b) {
        assert(isInSameSubsystem(b));
        invalidateTopologyCache();

        const ConstrainedBodyId nextId((int)myConstrainedBodies.size());

        // Add to the Mobilized->Constrained map and check for duplicates.
        std::pair<Mobilized2ConstrainedMap::iterator, bool> result;
        result = myMobilizedBodies.insert(
            Mobilized2ConstrainedMap::value_type(b.getMobilizedBodyId(), nextId));
        assert(result.second); // can only add a body once

        // This is a new constrained body -- add it to the Conrained->Mobilized map too.
        myConstrainedBodies.push_back(b.getMobilizedBodyId());
        return nextId;
    }

    void realizeTopology(int& nxtQErr, int& nxtUErr, int& nxtMult) const;

    // These are measured from and expressed in the ancestor (A) frame.
    //TODO: should precalculate in State, return reference
    Transform  getBodyTransform   (const State& s, ConstrainedBodyId B) const; // X_AB
    SpatialVec getBodyVelocity    (const State& s, ConstrainedBodyId B) const; // V_AB
    SpatialVec getBodyAcceleration(const State& s, ConstrainedBodyId B) const; // A_AB

    // Extract just the rotational quantities from the spatial quantities above.
    const Rotation& getBodyRotation           (const State& s, ConstrainedBodyId B) const {return getBodyTransform(s,B).R();}   // R_AB
    const Vec3&     getBodyAngularVelocity    (const State& s, ConstrainedBodyId B) const {return getBodyVelocity(s,B)[0];}     // w_AB
    const Vec3&     getBodyAngularAcceleration(const State& s, ConstrainedBodyId B) const {return getBodyAcceleration(s,B)[0];} // b_AB

    // Extract just the translational (linear) quantities from the spatial quantities above.
    const Vec3& getBodyOriginLocation    (const State& s, ConstrainedBodyId B) const {return getBodyTransform(s,B).T();}   // p_AB
    const Vec3& getBodyOriginVelocity    (const State& s, ConstrainedBodyId B) const {return getBodyVelocity(s,B)[1];}     // v_AB
    const Vec3& getBodyOriginAcceleration(const State& s, ConstrainedBodyId B) const {return getBodyAcceleration(s,B)[1];} // a_AB

    Vec3 calcStationLocation(const State& s, ConstrainedBodyId B, const Vec3& p_B) const {
        return getBodyTransform(s,B) * p_B; // re-measure and re-express
    }
    Vec3 calcStationVelocity(const State& s, ConstrainedBodyId B, const Vec3& p_B) const {
        const Vec3 p_A = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const SpatialVec& V_AB = getBodyVelocity(s,B);
        return V_AB[1] + V_AB[0] % p_A;
    }
    Vec3 calcStationAcceleration(const State& s, ConstrainedBodyId B, const Vec3& p_B) const {
        const Vec3 p_A = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const SpatialVec& V_AB = getBodyVelocity(s,B);
        const SpatialVec& A_AB = getBodyAcceleration(s,B);
        return A_AB[1] + A_AB[0] % p_A + V_AB[0] % (V_AB[1] + V_AB[0] % p_A);
    }

    virtual ~ConstraintRep();
    virtual ConstraintRep* clone() const = 0;

    // This creates a constraint node using the appropriate constraint type.
    virtual ConstraintNode* createConstraintNode() const = 0; 


    // After realizeTopology() we can look at the values of modeling variables in the State.
    // A Constraint is free to use those in determining how many constraint equations of each
    // type to generate. The default implementation here doesn't look at the state but instead
    // returns the default numbers of equations supplied when the Constraint was constructed.
    void calcNumConstraintEquations(const State& s, int& mp, int& mv, int& ma) const {
        calcNumConstraintEquationsVirtual(s,mp,mv,ma);
    }
    virtual void calcNumConstraintEquationsVirtual(const State&, int& mp, int& mv, int& ma) const {
        mp = defaultMp; mv = defaultMv; ma = defaultMa;
    }

    //NOTE: bodyForces and mobilityForces refer only to constrained bodies and their
    //associated mobilizers, not the system as a whole. They are initialized to zero
    //prior to the call so do not need to be set.

    //NOTE: each of these operators acts on the current state of this constraint's
    //Subtree, which may or may not be the same as that Subtree has in the global
    //State. This is controlled by the base class operator interface methods which
    //will call these only after setting the Subtree state properly.
    //TODO: Subtree

    void calcPositionErrors(const State& s, int mp,  Real* perr) const {
        calcPositionErrorsVirtual(s,mp,perr);
    }
    void calcPositionDotErrors(const State& s, int mp,  Real* pverr) const {
        calcPositionDotErrorsVirtual(s,mp,pverr);
    }
    void calcPositionDotDotErrors(const State& s, int mp,  Real* paerr) const {
        calcPositionDotDotErrorsVirtual(s,mp,paerr);
    }
    void applyPositionConstraintForces
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyPositionConstraintForcesVirtual(s,mp,multipliers,bodyForces,mobilityForces);
    }

    void calcVelocityErrors(const State& s, int mv,  Real* verr) const {
        calcVelocityErrorsVirtual(s,mv,verr);
    }
    void calcVelocityDotErrors(const State& s, int mv,  Real* vaerr) const {
        calcVelocityDotErrorsVirtual(s,mv,vaerr);
    }
    void applyVelocityConstraintForces
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyVelocityConstraintForcesVirtual(s,mv,multipliers,bodyForces,mobilityForces);
    }

    void calcAccelerationErrors(const State& s, int ma,  Real* aerr) const {
        calcAccelerationErrorsVirtual(s,ma,aerr);
    }
    void applyAccelerationConstraintForces
       (const State& s, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyAccelerationConstraintForcesVirtual(s,ma,multipliers,bodyForces,mobilityForces);
    }

    // These must be defined if there are any position (holonomic) constraints defined.
    virtual void calcPositionErrorsVirtual      (const State&, int mp,  Real* perr) const;
    virtual void calcPositionDotErrorsVirtual   (const State&, int mp,  Real* pverr) const;
    virtual void calcPositionDotDotErrorsVirtual(const State&, int mp,  Real* paerr) const;
    virtual void applyPositionConstraintForcesVirtual
       (const State&, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any velocity (nonholonomic) constraints defined.
    virtual void calcVelocityErrorsVirtual   (const State&, int mv,  Real* verr) const;
    virtual void calcVelocityDotErrorsVirtual(const State&, int mv,  Real* vaerr) const;
    virtual void applyVelocityConstraintForcesVirtual
       (const State&, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any acceleration-only constraints defined.
    virtual void calcAccelerationErrorsVirtual(const State&, int ma,  Real* aerr) const;
    virtual void applyAccelerationConstraintForcesVirtual
       (const State&, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;


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
    bool subsystemTopologyHasBeenRealized() const;
    const ConstraintNode& getMyConstraintNode() const;

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              ConstraintId id);

    const SimbodyMatterSubsystem& getMyMatterSubsystem() const;

    bool isInSubsystem() const {
        return myMatterSubsystemRep != 0;
    }

    bool isInSameSubsystem(const MobilizedBody& body) const {
        return isInSubsystem() && body.isInSubsystem() 
            && getMyMatterSubsystem().isSameSubsystem(body.getMatterSubsystem());
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

    // We'll keep the bodies in two maps: one maps MobilizedBodyId->ConstrainedBodyId
    // (O(log n) to look up), and the other maps ConstrainedBodyId->MobilizedBodyId
    // (randomly addressable in constant time).
    Mobilized2ConstrainedMap     myMobilizedBodies;
    std::vector<MobilizedBodyId> myConstrainedBodies; // index with ConstrainedBodyId

    // These are the defaults for the number of position (holonomic) constraint equations,
    // the number of velocity (nonholonomic) constraint equations, and the number of
    // acceleration-only constraint equations.
    int defaultMp, defaultMv, defaultMa;

        // TOPOLOGY "CACHE"

    // When topology is realized we study the constrained bodies to identify the
    // subtree of mobilized bodies which may be kinematically involved in satisfaction
    // of this Constraint. This requires finding the outmost common ancestor of 
    // the constrained bodies. All mobilized bodies on the paths inward from the
    // constrained bodies to the ancestor are included; nothing outboard of the
    // constrained bodies is included; and the ancestor is treated as ground so
    // that its mobilities are *not* included.
    mutable SimbodyMatterSubsystem::Subtree mySubtree;

    // A ConstraintNode is created for each Constraint during realizeTopology().
    // Think of it as the *computational* form of the Constraint; the Constraint
    // itself exists to make a nice API for the programmer. In fact is is possible
    // for several different Constraint classes to use the same underlying computational
    // form while providing a different API.
    //
    // This is a pointer to an abstract object whose heap space is *owned* by
    // this Constraint. Be sure to delete it upon destruction and whenever
    // topology is re-realized.
    mutable ConstraintNode* myConstraintNode;
};


class Constraint::Rod::RodRep : public Constraint::ConstraintRep {
public:
    RodRep() 
      : ConstraintRep(1,0,0), defaultPoint1(0), defaultPoint2(0), defaultRodLength(1)
    { 
        // Rod constructor sets all the data members here directly
    }
    RodRep* clone() const { return new RodRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    // Implementation of virtuals required for holonomic constraints.

    // perr = (p^2 - d^2)/2
    void calcPositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==1 && perr);
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        //TODO: save p in state

        *perr = (dot(p, p) - square(defaultRodLength)) / 2;
    }

    // pverr = d/dt perr = pdot*p = v*p, where v=v2-v1 is relative velocity
    void calcPositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p from State
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;

        const Vec3 v1 = calcStationVelocity(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = calcStationVelocity(s, B2, defaultPoint2);
        const Vec3 v = v2 - v1;
        *pverr = dot(v, p);
    }

    // paerr = d/dt verr = vdot*p + v*pdot =a*p+v*v, where a=a2-a1 is relative acceleration
    void calcPositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v from State
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        const Vec3 v1 = calcStationVelocity(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = calcStationVelocity(s, B2, defaultPoint2);
        const Vec3 v = v2 - v1;

        const Vec3 a1 = calcStationAcceleration(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 a2 = calcStationAcceleration(s, B2, defaultPoint2);
        const Vec3 a = a2 - a1;

        *paerr = dot(a, p) + dot(v, v);
    }

    // Write this routine by inspection of the pverr routine, looking for terms involving
    // velocity. On point2 we see v2*p, on point1 we see -v1*p, so forces are m*p and -m*p,
    // respectively.
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        assert(mp==1 && multipliers);
        const Real lambda = *multipliers;
        //TODO: should be able to get p from State
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;

        bodyForces[B2][1] = lambda * p; // no torque
        bodyForces[B1][1] = -bodyForces[B2][1];
    }

    SimTK_DOWNCAST(RodRep, ConstraintRep);
private:
    friend class Constraint::Rod;

    ConstrainedBodyId B1, B2;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRodLength;
};


class Constraint::PointInPlane::PointInPlaneRep : public Constraint::ConstraintRep {
public:
    PointInPlaneRep()
      : ConstraintRep(1,0,0), defaultPlaneNormal(), defaultPlaneHeight(0), defaultFollowerPoint(0),
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
    BallRep() : ConstraintRep(3,0,0), defaultPoint1(0), defaultPoint2(0), defaultRadius(0.1) { }
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
    WeldRep() : ConstraintRep(6,0,0) { } // transforms are Identity
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



