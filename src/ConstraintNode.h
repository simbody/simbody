#ifndef SimTK_SIMBODY_CONSTRAINT_NODE_H_
#define SimTK_SIMBODY_CONSTRAINT_NODE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

#include "simbody/internal/common.h"
#include "SimbodyTreeState.h"

#include "LengthConstraints.h"

#include <cassert>
#include <vector>

using namespace SimTK;

/**
 * This class represents a "Constraint", which is in general a set of related
 * constraint equations. 
 */
class ConstraintNode {
public:
    virtual ~ConstraintNode() {}

    ConstraintNode& operator=(const ConstraintNode&);

    virtual const char* type()      const=0;
    virtual ConstraintNode* clone() const=0;


        // TOPOLOGICAL INFO: no State needed

    int  getConstraintNum() const {return constraintNum;}
    void setConstraintNum(int n)  {constraintNum=n;}

    int  getQErrIndex() const {assert(qerrIndex>=0); return qerrIndex;}
    void setQErrIndex(int ix) {qerrIndex=ix;}

    int  getUErrIndex() const {assert(uerrIndex>=0); return uerrIndex;}
    void setUErrIndex(int ix) {uerrIndex=ix;}

    int  getMultIndex() const {assert(multIndex>=0); return multIndex;}
    void setMultIndex(int ix) {multIndex=ix;}

    // TODO: should this be variable based on modeling stage info?
    virtual int getNConstraintEquations() const=0; // # lambda slots

    
    // This gives each constraint node a chance to grab resources from
    // the SimbodyMatterSubsystemRep. The node is responsible for remembering which
    // resources belong to it. Currently the only resource is a set
    // of distance constraints and the corresponding state entries.
    // We expect to know our allocated multipliers prior to this call so
    // we can dole out the slots to our subsidiary constraint equations.
    virtual void finishConstruction(SimbodyMatterSubsystemRep&) {assert(false);}


protected:
    // This is the constructor for the abstract base type for use by the derived
    // concrete types in their constructors.
    ConstraintNode() : qerrIndex(-1), uerrIndex(-1), multIndex(-1), constraintNum(-1) { }

    int qerrIndex, uerrIndex, multIndex;
    int constraintNum;  // unique ID number in SimbodyMatterSubsystemRep

    friend std::ostream& operator<<(std::ostream& s, const ConstraintNode&);
};

/**
 * This class represents a single constraint equation, maintaining a fixed distance
 * between stations on two different bodies (or a body and ground).
 */
class ConstantDistanceConstraintNode : public ConstraintNode {
public:
    ConstantDistanceConstraintNode(
            const RigidBodyNode& parent, const Vec3& stationInP,
            const RigidBodyNode& child,  const Vec3& stationInC,
            const Real& distance)
      : body1(parent), body2(child), station1(stationInP), station2(stationInC),
        separation(distance), distanceConstraintIndex(-1)
    {
    }
    ~ConstantDistanceConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
        const RBStation s1(body1, station1);
        const RBStation s2(body2, station2);
        distanceConstraintIndex = 
            tree.addOneDistanceConstraintEquation(s1,s2,separation,
                    getQErrIndex(), getUErrIndex(), getMultIndex());
    }

    /*virtual*/ const char* type()     const {return "constantDistance";}
    /*virtual*/ int         getNConstraintEquations() const {return 1;}
    /*virtual*/ ConstraintNode* clone() const {
        return new ConstantDistanceConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& body1;
    const RigidBodyNode& body2;  
    const Vec3 station1, station2;
    const Real separation;

    int distanceConstraintIndex;
};

/**
 * This class represents a single constraint equation, enforcing that a point fixed
 * to one body move in a plane fixed to another.
 */
class PointInPlaneConstraintNode : public ConstraintNode {
public:
    PointInPlaneConstraintNode(
            const RigidBodyNode& planeNode, const UnitVec3& planeNormalInP, Real planeHeight,
            const RigidBodyNode& followerNode,  const Vec3& followerPointInF)
      : planeBody(planeNode), followerBody(followerNode), normal(planeNormalInP), height(planeHeight),
        followerStation(followerPointInF), pointInPlaneConstraintIndex(-1)
    {
    }
    ~PointInPlaneConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
        const RBDirection d(planeBody, normal);
        const RBStation   s(followerBody, followerStation);
        pointInPlaneConstraintIndex = 
            tree.addOnePointInPlaneEquation(d, height, s,
                    getQErrIndex(), getUErrIndex(), getMultIndex());
    }

    /*virtual*/ const char* type()     const {return "pointInPlane";}
    /*virtual*/ int         getNConstraintEquations() const {return 1;}
    /*virtual*/ ConstraintNode* clone() const {
        return new PointInPlaneConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& planeBody;
    const RigidBodyNode& followerBody;  
    const UnitVec3 normal;
    const Real     height;
    const Vec3     followerStation;

    int pointInPlaneConstraintIndex;
};

/**
 * This class represents two constraint equations, enforcing that a point fixed
 * to one body move on a line fixed to another.
 * TODO: THIS IS A USELESS STUB TO BE REMOVED
 */
class PointOnLineConstraintNode : public ConstraintNode {
public:
    PointOnLineConstraintNode(
            const RigidBodyNode& lineNode, const UnitVec3& lineDirectionInB, const Vec3& pointOnLineInB,
            const RigidBodyNode& followerNode,  const Vec3& followerPointInF)
      : lineBody(lineNode), followerBody(followerNode), dir(lineDirectionInB), point(pointOnLineInB),
        followerStation(followerPointInF), pointOnLineConstraintIndex(-1)
    {
    }
    ~PointOnLineConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
    }

    /*virtual*/ const char* type()     const {return "pointOnLine";}
    /*virtual*/ int         getNConstraintEquations() const {return 2;}
    /*virtual*/ ConstraintNode* clone() const {
        return new PointOnLineConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& lineBody;
    const RigidBodyNode& followerBody;  
    const UnitVec3 dir;
    const Vec3     point;
    const Vec3     followerStation;

    int pointOnLineConstraintIndex;
};

/**
 * This class represents a single constraint equation, enforcing that a vector fixed
 * to one body maintains a constant angle with respect to a vector fixed to another
 * body.
 */
class ConstantAngleConstraintNode : public ConstraintNode {
public:
    ConstantAngleConstraintNode(
            const RigidBodyNode& baseNode, const UnitVec3& axisOnB,
            const RigidBodyNode& followerNode,  const UnitVec3& axisOnF, Real a)
      : baseBody(baseNode), followerBody(followerNode), 
      axisB(axisOnB), axisF(axisOnF), angle(a), constantAngleConstraintIndex(-1)
    {
    }
    ~ConstantAngleConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
    }

    /*virtual*/ const char* type()     const {return "constantAngle";}
    /*virtual*/ int         getNConstraintEquations() const {return 1;}
    /*virtual*/ ConstraintNode* clone() const {
        return new ConstantAngleConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& baseBody;
    const RigidBodyNode& followerBody;  
    const UnitVec3 axisB;
    const UnitVec3 axisF;
    const Real     angle;

    int constantAngleConstraintIndex;
};

/**
 * This class represents three constraint equations, conspiring to keep a frame on
 * one body aligned with one on another body.
 */
class ConstantOrientationConstraintNode : public ConstraintNode {
public:
    ConstantOrientationConstraintNode(
            const RigidBodyNode& baseNode, const Rotation& frameOnB,
            const RigidBodyNode& followerNode,  const Rotation& frameOnF)
      : baseBody(baseNode), followerBody(followerNode), 
      frameB(frameOnB), frameF(frameOnF), constantOrientationConstraintIndex(-1)
    {
    }
    ~ConstantOrientationConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
    }

    /*virtual*/ const char* type()     const {return "constantOrientation";}
    /*virtual*/ int         getNConstraintEquations() const {return 3;}
    /*virtual*/ ConstraintNode* clone() const {
        return new ConstantOrientationConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& baseBody;
    const RigidBodyNode& followerBody;  
    const Rotation frameB;
    const Rotation frameF;

    int constantOrientationConstraintIndex;
};

/**
 * This class represents three constraint equations, together conspiring to hold
 * a station on one body coincident with one on another (or on ground).
 * The current implementation uses three distance constraints between the
 * pair of bodies. TODO: that doesn't work very well! (0,0,0) and (2/3,2/3,2/3)
 * both satisfy the constraint equations generated below (each is exactly 1
 * unit away from all points 100, 010, 001).
 */
class CoincidentStationsConstraintNode : public ConstraintNode {
public:
    CoincidentStationsConstraintNode(
            const RigidBodyNode& parent, const Vec3& stationInP,
            const RigidBodyNode& child,  const Vec3& stationInC)
      : body1(parent), body2(child), station1(stationInP), station2(stationInC),
        firstDistanceConstraintIndex(-1)
    {
    }
    ~CoincidentStationsConstraintNode() { }

    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
        const RBStation s1x(body1, station1+Vec3(1,0,0));
        const RBStation s1y(body1, station1+Vec3(0,1,0));       
        const RBStation s1z(body1, station1+Vec3(0,0,1));
        const RBStation s2(body2, station2);
        const int qx = getQErrIndex();
        const int ux = getUErrIndex();
        const int mx = getMultIndex();
        firstDistanceConstraintIndex = 
              tree.addOneDistanceConstraintEquation(s1x,s2,1.,qx+0,ux+0,mx+0);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2,1.,qx+1,ux+1,mx+1);
        (void)tree.addOneDistanceConstraintEquation(s1z,s2,1.,qx+2,ux+2,mx+2);
    }

    /*virtual*/ const char* type()     const {return "separation";}
    /*virtual*/ int         getNConstraintEquations() const {return 3;}
    /*virtual*/ ConstraintNode* clone() const {
        return new CoincidentStationsConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& body1;
    const RigidBodyNode& body2;  
    const Vec3 station1, station2;

    int firstDistanceConstraintIndex;
};

/**
 * This class represents six constraint equations, together serving to weld
 * one body to another (or to ground).
 */
class WeldConstraintNode : public ConstraintNode {
public:
    WeldConstraintNode(
            const RigidBodyNode& parent, const Transform& frameInP,
            const RigidBodyNode& child,  const Transform& frameInC)
      : body1(parent), body2(child), frame1(frameInP), frame2(frameInC),
        firstDistanceConstraintIndex(-1)
    {
    }
    ~WeldConstraintNode() { }


    /*virtual*/ void finishConstruction(SimbodyMatterSubsystemRep& tree) {
        const Vec3& station1 = frame1.T();
        const Vec3& station2 = frame2.T();

        const RBStation s1 (body1, station1);
        const RBStation s1x(body1, station1+frame1.x());
        const RBStation s1y(body1, station1+frame1.y());       
        const RBStation s1z(body1, station1+frame1.z());

        const RBStation s2 (body2, station2);
        const RBStation s2x(body2, station2+frame2.x());
        const RBStation s2y(body2, station2+frame2.y());       
        const RBStation s2z(body2, station2+frame2.z());

        const int qx = getQErrIndex();
        const int ux = getUErrIndex();
        const int mx = getMultIndex();

        // This is a "coincident station" constraint holding the frame
        // origins together (see above).
        firstDistanceConstraintIndex = 
              tree.addOneDistanceConstraintEquation(s1x,s2,1.,qx+0,ux+0,mx+0);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2,1.,qx+1,ux+1,mx+1);
        (void)tree.addOneDistanceConstraintEquation(s1z,s2,1.,qx+2,ux+2,mx+2);

        // This is an "align axes" constraint. This has an unfortunate
        // symmetry when rotating 180 degrees about any axis.
        // This set of constraint equations is fine for *projection* but
        // not enough for *assembly*. You need to add another one to
        // eliminate the rotational symmetries when assembling from
        // far away.
        const Real d = std::sqrt(2.);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2z,d,qx+3,ux+3,mx+3); // restrain x rot
        (void)tree.addOneDistanceConstraintEquation(s1z,s2x,d,qx+4,ux+4,mx+4); // restrain y rot
        (void)tree.addOneDistanceConstraintEquation(s1x,s2y,d,qx+5,ux+5,mx+5); // restrain z rot    
    }

    /*virtual*/ const char* type()     const {return "separation";}
    /*virtual*/ int         getNConstraintEquations() const {return 6;}
    /*virtual*/ ConstraintNode* clone() const {
        return new WeldConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& body1;
    const RigidBodyNode& body2;  
    const Transform frame1, frame2;

    int firstDistanceConstraintIndex;
};

#endif // SimTK_SIMBODY_CONSTRAINT_NODE_H_
