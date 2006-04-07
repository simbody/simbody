#ifndef SimTK_SIMBODY_CONSTRAINT_NODE_H_
#define SimTK_SIMBODY_CONSTRAINT_NODE_H_

#include "simbody/internal/common.h"
#include "SimbodyTreeState.h"

#include "LengthConstraints.h"

#include <cassert>
#include <vector>

using namespace SimTK;

/**
 * This class represents a "constraint", which is in general a set of related
 * constraint equations. 
 */
class ConstraintNode {
public:
    virtual ~ConstraintNode() {}

    ConstraintNode& operator=(const ConstraintNode&);

    virtual const char* type()     const=0;
    virtual ConstraintNode* clone() const=0;


        // TOPOLOGICAL INFO: no State needed

    int  getConstraintNum() const {return constraintNum;}
    void setConstraintNum(int n)  {constraintNum=n;}

    int  getMultIndex() const {assert(multIndex>=0); return multIndex;}
    void setMultIndex(int ix) {multIndex=ix;}

    // TODO: should this be variable based on modeling stage info?
    virtual int getNMult() const=0; // # lambda slots

    
    // This gives each constraint node a chance to grab resources from
    // the RigidBodyTree. The node is responsible for remembering which
    // resources belong to it. Currently the only resource is a set
    // of distance constraints and the corresponding state entries.
    // We expect to know our allocated multipliers prior to this call so
    // we can dole out the slots to our subsidiary constraint equations.
    virtual void finishConstruction(RigidBodyTree&) {assert(false);}


protected:
    // This is the constructor for the abstract base type for use by the derived
    // concrete types in their constructors.
    ConstraintNode() : multIndex(-1), constraintNum(-1) { }

    int               multIndex;      // index into lambda array
    int               constraintNum;  // unique ID number in RigidBodyTree

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

    /*virtual*/ void finishConstruction(RigidBodyTree& tree) {
        const RBStation s1(body1, station1);
        const RBStation s2(body2, station2);
        distanceConstraintIndex = 
            tree.addOneDistanceConstraintEquation(s1,s2,separation,getMultIndex());
    }

    /*virtual*/ const char* type()     const {return "separation";}
    /*virtual*/ int         getNMult() const {return 1;}
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
 * This class represents three constraint equations, together conspiring to hold
 * a station on one body coincident with one on another (or on ground).
 * The current implementation uses three distance constraints between the
 * pair of bodies.
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

    /*virtual*/ void finishConstruction(RigidBodyTree& tree) {
        const RBStation s1x(body1, station1+Vec3(1,0,0));
        const RBStation s1y(body1, station1+Vec3(0,1,0));       
        const RBStation s1z(body1, station1+Vec3(0,0,1));
        const RBStation s2(body2, station2);
        const int ix = getMultIndex();
        firstDistanceConstraintIndex = 
              tree.addOneDistanceConstraintEquation(s1x,s2,1.,ix+0);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2,1.,ix+1);
        (void)tree.addOneDistanceConstraintEquation(s1z,s2,1.,ix+2);
    }

    /*virtual*/ const char* type()     const {return "separation";}
    /*virtual*/ int         getNMult() const {return 3;}
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
            const RigidBodyNode& parent, const TransformMat& frameInP,
            const RigidBodyNode& child,  const TransformMat& frameInC)
      : body1(parent), body2(child), frame1(frameInP), frame2(frameInC),
        firstDistanceConstraintIndex(-1)
    {
    }
    ~WeldConstraintNode() { }


    /*virtual*/ void finishConstruction(RigidBodyTree& tree) {
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

        const int ix = getMultIndex();

        // This is a "coincident station" constraint holding the frame
        // origins together (see above).
        firstDistanceConstraintIndex = 
              tree.addOneDistanceConstraintEquation(s1x,s2,1.,ix+0);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2,1.,ix+1);
        (void)tree.addOneDistanceConstraintEquation(s1z,s2,1.,ix+2);

        // This is an "align axes" constraint. This has an unfortunate
        // symmetry when rotating 180 degrees about any axis.
        // This set of constraint equations is fine for *projection* but
        // not enough for *assembly*. You need to add another one to
        // eliminate the rotational symmetries when assembling from
        // far away.
        const Real d = std::sqrt(2.);
        (void)tree.addOneDistanceConstraintEquation(s1y,s2z,d,ix+3); // restrain x rot
        (void)tree.addOneDistanceConstraintEquation(s1z,s2x,d,ix+4); // restrain y rot
        (void)tree.addOneDistanceConstraintEquation(s1x,s2y,d,ix+5); // restrain z rot    
    }

    /*virtual*/ const char* type()     const {return "separation";}
    /*virtual*/ int         getNMult() const {return 6;}
    /*virtual*/ ConstraintNode* clone() const {
        return new WeldConstraintNode(*this);
    }
    
private:
    const RigidBodyNode& body1;
    const RigidBodyNode& body2;  
    const TransformMat frame1, frame2;

    int firstDistanceConstraintIndex;
};

#endif // SimTK_SIMBODY_CONSTRAINT_NODE_H_
