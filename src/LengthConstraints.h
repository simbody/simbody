#ifndef SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_
#define SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_

#include "simbody/internal/common.h"
using namespace SimTK;

#include "RigidBodyTree.h"
#include "SimbodyTreeState.h"

#include "newtonRaphson.h"

#include <vector>
#include <iostream>
using std::ostream;

/**
 * A station is a point located on a particular rigid body. A station is
 * measured from the body frame origin and expressed in the body frame.
 */
class RBStation {
public:
    RBStation() : rbNode(0) { } // so we can have arrays of these
    RBStation(const RigidBodyNode& n, const Vec3& pos) 
        : rbNode(&n), station_B(pos) { }
    // default copy, assignment, destructor

    const RigidBodyNode& getNode()    const { assert(isValid()); return *rbNode; }
    const Vec3&          getPoint()   const { assert(isValid()); return station_B; }
    bool                 isValid()    const { return rbNode != 0; }
private:
    const RigidBodyNode* rbNode;
    Vec3                 station_B;
};
std::ostream& operator<<(std::ostream&, const RBStation&);

/**
 * This class requests that two stations, one on each of two rigid bodies,
 * be maintained at a certain separation distance at all times. This is 
 * an internal service provided by RigidBodyTrees and not something
 * built directly by users. User-requested Constraints may allocate one
 * or more distance constraints in the performance of their duties.
 *
 * Each distance constraint adds one constraint equation and is thus
 * associated with a particular constraint multiplier and matching
 * acceleration constraint error. However, not all multipliers come
 * from distance constraints, so we store both the distance constraint
 * number here and the multiplier index, for use in accessing the 
 * appropriate quantities in the state cache.
 */
class RBDistanceConstraint {
public:
    RBDistanceConstraint() : distance(-1.), distConstNum(-1), multIndex(-1) {}
    RBDistanceConstraint(const RBStation& s1, const RBStation& s2, const Real& d)
    {
        assert(s1.isValid() && s2.isValid() && d >= 0.);
        stations[0] = s1; stations[1] = s2; distance = d;
        distConstNum = multIndex = -1;
    }

    void calcPosInfo(
        SBConfigurationCache&       cc) const;
    void calcVelInfo(
        const SBConfigurationCache& cc, 
        SBMotionCache&              mc) const;
    void calcAccInfo(
        const SBConfigurationCache& cc, 
        const SBMotionCache&        mc,
        SBReactionCache&            rc) const;

    void setDistanceConstraintNum(int ix) {assert(ix>=0); distConstNum=ix;}
    int  getDistanceConstraintNum() const {assert(isValid()&&distConstNum>=0); return distConstNum;}

    void setMultIndex(int ix) {assert(ix>=0); multIndex=ix;}
    int  getMultIndex() const {assert(isValid()&&multIndex>=0); return multIndex;}

    // Currently these are not parameterizable so they have no state references.
    const Real&          getDistance()     const {return distance;}
    const RBStation&     getStation(int i) const {assert(isValid() && (i==1||i==2)); return stations[i-1];}
    const RigidBodyNode& getNode(int i)    const {return getStation(i).getNode();}
    const Vec3&          getPoint(int i)   const {return getStation(i).getPoint();}
    bool                 isValid()         const {return distance >= 0.;}

    // State access routines

        // CONFIGURED STAGE
    const Real& getPosErr(const SBConfigurationCache& cc) const {
        return cc.positionConstraintErrors[multIndex];
    }
    Real& updPosErr(SBConfigurationCache& cc) const {
        return cc.positionConstraintErrors[multIndex];
    }

    const Vec3& getStation_G(const SBConfigurationCache& cc, int i) const {
        assert(1 <= i && i <= 2);
        return cc.station_G[i-1][distConstNum];
    }
    Vec3& updStation_G(SBConfigurationCache& cc, int i) const {
        assert(1 <= i && i <= 2);
        return cc.station_G[i-1][distConstNum];
    }

    const Vec3& getPos_G(const SBConfigurationCache& cc, int i) const {
        assert(1 <= i && i <= 2);
        return cc.pos_G[i-1][distConstNum];
    }
    Vec3& updPos_G(SBConfigurationCache& cc, int i) const {
        assert(1 <= i && i <= 2);
        return cc.pos_G[i-1][distConstNum];
    }

    const Vec3& getFromTip1ToTip2_G(const SBConfigurationCache& cc) const {
        return cc.fromTip1ToTip2_G[distConstNum];
    }
    Vec3& updFromTip1ToTip2_G(SBConfigurationCache& cc) const {
        return cc.fromTip1ToTip2_G[distConstNum];
    }

    const Vec3& getUnitDirection_G(const SBConfigurationCache& cc) const {
        return cc.unitDirection_G[distConstNum];
    }
    Vec3& updUnitDirection_G(SBConfigurationCache& cc) const {
        return cc.unitDirection_G[distConstNum];
    }

        // MOVING STAGE
    const Real& getVelErr(const SBMotionCache& mc) const {
        return mc.velocityConstraintErrors[multIndex];
    }
    Real& updVelErr(SBMotionCache& mc) const {
        return mc.velocityConstraintErrors[multIndex];
    }

    const Vec3& getStationVel_G(const SBMotionCache& mc, int i) const {
        assert(1 <= i && i <= 2);
        return mc.stationVel_G[i-1][distConstNum];
    }
    Vec3& updStationVel_G(SBMotionCache& mc, int i) const {
        assert(1 <= i && i <= 2);
        return mc.stationVel_G[i-1][distConstNum];
    }

    const Vec3& getVel_G(const SBMotionCache& mc, int i) const {
        assert(1 <= i && i <= 2);
        return mc.vel_G[i-1][distConstNum];
    }
    Vec3& updVel_G(SBMotionCache& mc, int i) const {
        assert(1 <= i && i <= 2);
        return mc.vel_G[i-1][distConstNum];
    }

    const Vec3& getRelVel_G(const SBMotionCache& mc) const {
        return mc.relVel_G[distConstNum];
    }
    Vec3& updRelVel_G(SBMotionCache& mc) const {
        return mc.relVel_G[distConstNum];
    }


        // REACTING STAGE
    const Real& getAccErr(const SBReactionCache& rc) const {
        return rc.accelerationConstraintErrors[multIndex];
    }
    Real& updAccErr(SBReactionCache& rc) const {
        return rc.accelerationConstraintErrors[multIndex];
    }

    const Vec3& getAcc_G(const SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        return rc.acc_G[i-1][distConstNum];
    }
    Vec3& updAcc_G(SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        return rc.acc_G[i-1][distConstNum];
    }

    const Vec3& getForce_G(const SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        return rc.force_G[i-1][distConstNum];
    }
    Vec3& updForce_G(SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        return rc.force_G[i-1][distConstNum];
    }

protected:
    Real       distance;
    RBStation  stations[2];
    int        distConstNum;
    int        multIndex;

private:
    // Per-station calculations
    void calcStationPosInfo(int i, 
        SBConfigurationCache&       cc) const;
    void calcStationVelInfo(int i, 
        const SBConfigurationCache& cc, 
        SBMotionCache&              mc) const;
    void calcStationAccInfo(int i, 
        const SBConfigurationCache& cc, 
        const SBMotionCache&        mc,
        SBReactionCache&            rc) const;
};


class RigidBodyTree;
class LengthConstraints;

/*
 * Collect up useful information about a loop. 
 * This includes the two connected stations, ordered by level, and the
 * paths from each of the associated nodes back to the common ancestor.
 * We also identify the molecule base node for the molecule which
 * contains both ends of the loop.
 * We will throw an exception if the loop ends are both on the same
 * node or if they are on different molecules.
 */
class LoopWNodes {
public:
    LoopWNodes() 
      : tree(0), rbDistCons(0), flipStations(false), outmostCommonBody(0) 
    {
    }
    LoopWNodes(const RigidBodyTree&, const RBDistanceConstraint&);

    void calcPosInfo(
        SBConfigurationCache&       cc) const 
      { rbDistCons->calcPosInfo(cc); }

    void calcVelInfo(
        const SBConfigurationCache& cc, 
        SBMotionCache&              mc) const 
      { rbDistCons->calcVelInfo(cc,mc); }

    void calcAccInfo(
        const SBConfigurationCache& cc, 
        const SBMotionCache&        mc,
        SBReactionCache&            rc) const 
      { rbDistCons->calcAccInfo(cc,mc,rc); }

    // TODO: State isn't currently in use but will be necessary when the distances
    // are parameters.
    const Real& getDistance() const { return rbDistCons->getDistance(); }

    // Return one of the stations, ordered such that tips(1).level <= tips(2).level.
    const RBStation& tips(int i) const {return rbDistCons->getStation(ix(i));}
    const RigidBodyNode& tipNode(int i) const {return tips(i).getNode();}
 
    const Vec3& tipPos(const SBConfigurationCache& cc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return cc.pos_G[ix(i)-1][dc];
    }
    const Vec3& tipVel(const SBMotionCache& mc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return mc.vel_G[ix(i)-1][dc];
    }
    const Vec3& tipAcc(const SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return rc.acc_G[ix(i)-1][dc];
    }
    const Vec3& tipForce(const SBReactionCache& rc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return rc.force_G[ix(i)-1][dc];
    }

    // Use this for both forces and impulses.
    void setTipForce(SBReactionCache& rc, int i, const Vec3& f) const {
        assert(1 <= i && i <= 2);

        const int dc = rbDistCons->getDistanceConstraintNum();
        rc.force_G[ix(i)-1][dc] = f;
    }

    const RigidBodyNode* getOutmostCommonBody() const {return outmostCommonBody;}

private:
    int ix(int i) const { assert(i==1||i==2); return flipStations ? 3-i : i; }

    const RigidBodyTree*         tree;        // a reference to the tree we're part of
    const RBDistanceConstraint*  rbDistCons;  // reference to the constraint

    // calculated construction-time (topological) info about the constraint
    bool                              flipStations; // make sure station(1).level
                                                    //   <= station(2).level
    std::vector<const RigidBodyNode*> nodes[2];     // the two paths: base..tip1, base..tip2,
                                                    //   incl. tip nodes but not base
    const RigidBodyNode*              outmostCommonBody; // highest-level common ancestor of tips
    
    // Ancestors includes everything from outmostCommonBody (inclusive) down to ground
    // (exclusive).
    std::vector<const RigidBodyNode*> ancestors;

    friend class LengthSet;
    friend class LengthConstraints;
    friend ostream& operator<<(ostream&, const LengthSet&);
    friend ostream& operator<<(ostream&, const LoopWNodes&);
};

typedef std::vector<LoopWNodes> LoopList;

class LengthSet {
    static void construct(const LoopList& loops);
    const LengthConstraints*          lConstraints;
    LoopList                          loops;    
    int                               ndofThisSet;
    std::vector<const RigidBodyNode*> nodeMap; //unique nodes (union of loops->nodes)
public:
    LengthSet() : lConstraints(0), ndofThisSet(0) { }
    LengthSet(const LengthConstraints* lConstraints)
      : lConstraints(lConstraints), ndofThisSet(0)
    {
    }

    inline const RigidBodyTree& getRBTree()  const;
    inline int                  getVerbose() const;

    void addKinematicConstraint(const LoopWNodes& loop);
    void addDynamicConstraint(const LoopWNodes& loop);
    bool contains(const RigidBodyNode* node);

    void  setPos(State&, const Vector& pos) const;
    void  setVel(State&, const Vector& vel) const;
    Vector getPos();
    Vector calcPosB(State&, const Vector& pos) const;
    Vector calcVelB(State&, const Vector& vel) const;
    Vector calcPosZ(const State&, const Vector& b) const;
    Matrix calcGrad(const State&) const;
    Matrix calcGInverse(const State&) const;
    Matrix calcGInverseFD(const State&) const;

    void  calcConstraintForces(const State&) const; // updates runtime only
    void  addInCorrectionForces(const State&, SpatialVecList& spatialForces) const; // spatialForces+=correction

    void   fixVel0(State&, Vector&);
    void   fixInternalForce(const State&, Vector&);
    Vector multiForce(const Vector&, const Matrix& mat);
    void   subtractVecFromForces(Vector& frc, const Vector& vec);

    void testAccel(const State&) const;
    void testInternalForce(const State&, const Vector&) const;
    friend ostream& operator<<(ostream& os, const LengthSet& s);
    friend class CalcPosB;
    friend class CalcVelB;
    friend class CalcPosZ;
    friend class CalcVelZ;

    void fdgradf(State&, const Vector& pos, Matrix& grad) const;
    void testGrad(State&, const Vector& pos, const Matrix& grad) const;
};

class LengthConstraints {
public:
    LengthConstraints(const RigidBodyTree&, int verbose);

    void construct(const Array<RBDistanceConstraint*>&);

    // Returns true if any change was made in the state.
    bool enforceConfigurationConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;
    bool enforceMotionConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    bool calcConstraintForces(const State&) const;
    void addInCorrectionForces(const State&, SpatialVecList& spatialForces) const;

    void fixVel0(State&, Vector&);
    void fixGradient(const State&, Vector&);

private:
    int    maxIters;
    int    maxMin;

    const RigidBodyTree&       rbTree;
    const int                  verbose;

    std::vector<LengthSet> pvConstraints;   // used for pos, vel
    std::vector<LengthSet> accConstraints;  // used for acc
    NewtonRaphson          posMin, velMin;

    friend class LengthSet;
};

inline const RigidBodyTree& 
LengthSet::getRBTree()  const {return lConstraints->rbTree;}
inline int                  
LengthSet::getVerbose() const {return lConstraints->verbose;}


#endif // SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_


