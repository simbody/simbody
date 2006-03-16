#ifndef LENGTH_CONSTRAINTS_H_
#define LENGTH_CONSTRAINTS_H_

#include "simbody/Simbody.h"
using namespace simtk;

#include "RigidBodyTree.h"
#include "SimbodyTreeState.h"

#include <vector>

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
 * built directly by users. User-requested constraints may allocate one
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
    RBDistanceConstraint(const RBStation& s1, const RBStation& s2, const Real& d) {
        assert(s1.isValid() && s2.isValid() && d >= 0.);
        stations[0] = s1; stations[1] = s2; distance = d;
        distConstNum = multIndex = -1;
    }

    void calcPosInfo(const SBStateRep&) const;
    void calcVelInfo(const SBStateRep&) const;
    void calcAccInfo(const SBStateRep&) const;

    void setDistanceConstraintNum(int ix) {assert(ix>=0); distConstNum=ix;}
    int  getDistanceConstraintNum() const {assert(isValid()&&distConstNum>=0); return distConstNum;}

    void setMultIndex(int ix) {assert(ix>=0); multIndex=ix;}
    int  getMultIndex() const {assert(isValid()&&multIndex>=0); return multIndex;}

    const Real&          getDistance()     const {return distance;}
    const RBStation&     getStation(int i) const {assert(isValid() && (i==1||i==2)); return stations[i-1];}
    const RigidBodyNode& getNode(int i)    const {return getStation(i).getNode();}
    const Vec3&          getPoint(int i)   const {return getStation(i).getPoint();}
    bool                 isValid()         const {return distance >= 0.;}

    // State access routines

        // CONFIGURED STAGE
    const Real& getPosErr(const SBStateRep& s) const {
        return s.configCache.positionConstraintErrors[multIndex];
    }
    Real& updPosErr(const SBStateRep& s) const {
        return s.configCache.positionConstraintErrors[multIndex];
    }

    const Vec3& getStation_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.configCache.station_G[i-1][distConstNum];
    }
    Vec3& updStation_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.configCache.station_G[i-1][distConstNum];
    }

    const Vec3& getPos_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.configCache.pos_G[i-1][distConstNum];
    }
    Vec3& updPos_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.configCache.pos_G[i-1][distConstNum];
    }

    const Vec3& getFromTip1ToTip2_G(const SBStateRep& s) const {
        return s.configCache.fromTip1ToTip2_G[distConstNum];
    }
    Vec3& updFromTip1ToTip2_G(const SBStateRep& s) const {
        return s.configCache.fromTip1ToTip2_G[distConstNum];
    }

    const Vec3& getUnitDirection_G(const SBStateRep& s) const {
        return s.configCache.unitDirection_G[distConstNum];
    }
    Vec3& updUnitDirection_G(const SBStateRep& s) const {
        return s.configCache.unitDirection_G[distConstNum];
    }

        // MOVING STAGE
    const Real& getVelErr(const SBStateRep& s) const {
        return s.motionCache.velocityConstraintErrors[multIndex];
    }
    Real& updVelErr(const SBStateRep& s) const {
        return s.motionCache.velocityConstraintErrors[multIndex];
    }

    const Vec3& getStationVel_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.motionCache.stationVel_G[i-1][distConstNum];
    }
    Vec3& updStationVel_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.motionCache.stationVel_G[i-1][distConstNum];
    }

    const Vec3& getVel_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.motionCache.vel_G[i-1][distConstNum];
    }
    Vec3& updVel_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.motionCache.vel_G[i-1][distConstNum];
    }

    const Vec3& getRelVel_G(const SBStateRep& s) const {
        return s.motionCache.relVel_G[distConstNum];
    }
    Vec3& updRelVel_G(const SBStateRep& s) const {
        return s.motionCache.relVel_G[distConstNum];
    }


        // REACTING STAGE
    const Real& getAccErr(const SBStateRep& s) const {
        return s.reactionCache.accelerationConstraintErrors[multIndex];
    }
    Real& updAccErr(const SBStateRep& s) const {
        return s.reactionCache.accelerationConstraintErrors[multIndex];
    }

    const Vec3& getAcc_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.reactionCache.acc_G[i-1][distConstNum];
    }
    Vec3& updAcc_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.reactionCache.acc_G[i-1][distConstNum];
    }

    const Vec3& getForce_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.reactionCache.force_G[i-1][distConstNum];
    }
    Vec3& updForce_G(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        return s.reactionCache.force_G[i-1][distConstNum];
    }

protected:
    Real       distance;
    RBStation  stations[2];
    int        distConstNum;
    int        multIndex;

private:
    // Per-station calculations
    void calcStationPosInfo(const SBStateRep& s, int i) const;
    void calcStationVelInfo(const SBStateRep& s, int i) const;
    void calcStationAccInfo(const SBStateRep& s, int i) const;
};


class RigidBodyTree;

class LengthConstraintsPrivates;
class LengthSet;

class LengthConstraints {
public:
    LengthConstraints(const RigidBodyTree&, const double& ctol, int verbose);
    ~LengthConstraints();
    void construct(const std::vector<RBDistanceConstraint*>&);

    // Returns true if any change was made in the state.
    bool enforceConfigurationConstraints(SBStateRep&) const;
    bool enforceMotionConstraints(SBStateRep&) const;

    bool calcConstraintForces(const SBStateRep&) const;
    void addInCorrectionForces(const SBStateRep&, SpatialVecList& spatialForces) const;

    void fixVel0(SBStateRep&, Vector&);
    void fixGradient(const SBStateRep&, Vector&);

private:
    //  double tol;
    double bandCut;  //cutoff for calculation of constraint matrix
    int    maxIters;
    int    maxMin;

    const RigidBodyTree&       rbTree;
    const int                  verbose;
    LengthConstraintsPrivates* priv;

    friend class LengthSet;
};


#endif /* LENGTH_CONSTRAINTS_H_ */
