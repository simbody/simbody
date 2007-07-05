#ifndef SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_
#define SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
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

#include "simbody/internal/common.h"
using namespace SimTK;

#include "SimbodyMatterSubsystemRep.h"
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
 * an internal service provided by SimbodyMatterSubsystems and not something
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
    RBDistanceConstraint() : distance(-1.), distConstNum(-1), 
                             qerrIndex(-1), uerrIndex(-1), multIndex(-1) {}
    RBDistanceConstraint(const RBStation& s1, const RBStation& s2, const Real& d)
    {
        assert(s1.isValid() && s2.isValid() && d >= 0.);
        stations[0] = s1; stations[1] = s2; distance = d;
        distConstNum = qerrIndex = uerrIndex = multIndex = -1;
    }

    void calcPosInfo(
        Vector&                qErr,
        SBPositionCache&       pc) const;
    void calcVelInfo(
        const SBPositionCache& pc, 
        Vector&                uErr,
        SBVelocityCache&       vc) const;
    void calcAccInfo(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        Vector&                udotErr,
        SBAccelerationCache&   ac) const;

    void setDistanceConstraintNum(int ix) {assert(ix>=0); distConstNum=ix;}
    int  getDistanceConstraintNum() const {assert(isValid()&&distConstNum>=0); return distConstNum;}

    void setQErrIndex(int ix) {assert(ix>=0); qerrIndex=ix;}
    int  getQErrIndex() const {assert(isValid()&&qerrIndex>=0); return qerrIndex;}
    void setUErrIndex(int ix) {assert(ix>=0); uerrIndex=ix;}
    int  getUErrIndex() const {assert(isValid()&&uerrIndex>=0); return uerrIndex;}
    void setMultIndex(int ix) {assert(ix>=0); multIndex=ix;}
    int  getMultIndex() const {assert(isValid()&&multIndex>=0); return multIndex;}

    // Currently these are not parameterizable so they have no state references.
    const Real&          getDistance()     const {return distance;}
    const RBStation&     getStation(int i) const {assert(isValid() && (i==1||i==2)); return stations[i-1];}
    const RigidBodyNode& getNode(int i)    const {return getStation(i).getNode();}
    const Vec3&          getPoint(int i)   const {return getStation(i).getPoint();}
    bool                 isValid()         const {return distance >= 0.;}

    // State access routines

        // POSITION STAGE
    const Real& getPosErr(const Vector& qErr) const {
        assert(qerrIndex >= 0);
        return qErr[qerrIndex];
    }
    Real& updPosErr(Vector& qErr) const {
        assert(qerrIndex >= 0);
        return qErr[qerrIndex];
    }

    const Vec3& getStation_G(const SBPositionCache& pc, int i) const {
        assert(1 <= i && i <= 2);
        return pc.station_G[i-1][distConstNum];
    }
    Vec3& updStation_G(SBPositionCache& pc, int i) const {
        assert(1 <= i && i <= 2);
        return pc.station_G[i-1][distConstNum];
    }

    const Vec3& getPos_G(const SBPositionCache& pc, int i) const {
        assert(1 <= i && i <= 2);
        return pc.pos_G[i-1][distConstNum];
    }
    Vec3& updPos_G(SBPositionCache& pc, int i) const {
        assert(1 <= i && i <= 2);
        return pc.pos_G[i-1][distConstNum];
    }

    const Vec3& getFromTip1ToTip2_G(const SBPositionCache& pc) const {
        return pc.fromTip1ToTip2_G[distConstNum];
    }
    Vec3& updFromTip1ToTip2_G(SBPositionCache& pc) const {
        return pc.fromTip1ToTip2_G[distConstNum];
    }

    const Vec3& getUnitDirection_G(const SBPositionCache& pc) const {
        return pc.unitDirection_G[distConstNum];
    }
    Vec3& updUnitDirection_G(SBPositionCache& pc) const {
        return pc.unitDirection_G[distConstNum];
    }

        // VELOCITY STAGE
    const Real& getVelErr(const Vector& uErr) const {
        return uErr[uerrIndex];
    }
    Real& updVelErr(Vector& uErr) const {
        return uErr[uerrIndex];
    }

    const Vec3& getStationVel_G(const SBVelocityCache& vc, int i) const {
        assert(1 <= i && i <= 2);
        return vc.stationVel_G[i-1][distConstNum];
    }
    Vec3& updStationVel_G(SBVelocityCache& vc, int i) const {
        assert(1 <= i && i <= 2);
        return vc.stationVel_G[i-1][distConstNum];
    }

    const Vec3& getVel_G(const SBVelocityCache& vc, int i) const {
        assert(1 <= i && i <= 2);
        return vc.vel_G[i-1][distConstNum];
    }
    Vec3& updVel_G(SBVelocityCache& vc, int i) const {
        assert(1 <= i && i <= 2);
        return vc.vel_G[i-1][distConstNum];
    }

    const Vec3& getRelVel_G(const SBVelocityCache& vc) const {
        return vc.relVel_G[distConstNum];
    }
    Vec3& updRelVel_G(SBVelocityCache& vc) const {
        return vc.relVel_G[distConstNum];
    }


        // ACCELERATION STAGE
    const Real& getAccErr(const Vector& udotErr) const {
        return udotErr[multIndex];
    }
    Real& updAccErr(Vector& udotErr) const {
        return udotErr[multIndex];
    }

    const Vec3& getAcc_G(const SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        return ac.acc_G[i-1][distConstNum];
    }
    Vec3& updAcc_G(SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        return ac.acc_G[i-1][distConstNum];
    }

    const Vec3& getForce_G(const SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        return ac.force_G[i-1][distConstNum];
    }
    Vec3& updForce_G(SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        return ac.force_G[i-1][distConstNum];
    }

protected:
    Real       distance;
    RBStation  stations[2];
    int        distConstNum;
    int        qerrIndex, uerrIndex, multIndex;

private:
    // Per-station calculations
    void calcStationPosInfo(int i, 
        SBPositionCache&        pc) const;
    void calcStationVelInfo(int i, 
        const SBPositionCache&  pc, 
        SBVelocityCache&        vc) const;
    void calcStationAccInfo(int i, 
        const SBPositionCache&  pc, 
        const SBVelocityCache&  vc,
        SBAccelerationCache&    ac) const;
};


class SimbodyMatterSubsystemRep;
class LengthConstraints;
class LengthSet;

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
    LoopWNodes(const SimbodyMatterSubsystemRep&, const RBDistanceConstraint&);

    void calcPosInfo(
        Vector&          qErr,
        SBPositionCache& pc) const 
      { rbDistCons->calcPosInfo(qErr,pc); }

    void calcVelInfo(
        const SBPositionCache& pc, 
        Vector&                uErr,
        SBVelocityCache&       vc) const 
      { rbDistCons->calcVelInfo(pc,uErr,vc); }

    void calcAccInfo(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        Vector&                udotErr,
        SBAccelerationCache&   ac) const 
      { rbDistCons->calcAccInfo(pc,vc,udotErr,ac); }

    // TODO: State isn't currently in use but will be necessary when the distances
    // are parameters.
    const Real& getDistance() const { return rbDistCons->getDistance(); }

    // Return one of the stations, ordered such that tips(1).level <= tips(2).level.
    const RBStation& tips(int i) const {return rbDistCons->getStation(ix(i));}
    const RigidBodyNode& tipNode(int i) const {return tips(i).getNode();}
 
    const Vec3& tipPos(const SBPositionCache& pc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return pc.pos_G[ix(i)-1][dc];
    }
    const Vec3& tipVel(const SBVelocityCache& vc, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return vc.vel_G[ix(i)-1][dc];
    }
    const Vec3& tipAcc(const SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return ac.acc_G[ix(i)-1][dc];
    }
    const Vec3& tipForce(const SBAccelerationCache& ac, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return ac.force_G[ix(i)-1][dc];
    }

    // Use this for both forces and impulses.
    void setTipForce(SBAccelerationCache& ac, int i, const Vec3& f) const {
        assert(1 <= i && i <= 2);

        const int dc = rbDistCons->getDistanceConstraintNum();
        ac.force_G[ix(i)-1][dc] = f;
    }

    const RigidBodyNode* getOutmostCommonBody() const {return outmostCommonBody;}

private:
    int ix(int i) const { assert(i==1||i==2); return flipStations ? 3-i : i; }

    const SimbodyMatterSubsystemRep*         tree;        // a reference to the tree we're part of
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

    inline const SimbodyMatterSubsystemRep& getRBTree()  const;
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
    static Matrix calcPseudoInverseA(const Matrix& AT);
    Matrix calcPseudoInverseAFD(const State&) const;

    void  calcConstraintForces(const State&) const; // updates runtime only
    void  addInCorrectionForces(const State&, SpatialVecList& spatialForces) const; // spatialForces+=correction

    void   fixVel0(State&, Vector&);

    void   projectQVecOntoConfigurationConstraints(const State&, Vector& q);
    void   projectUVecOntoMotionConstraints(const State&, Vector& u);

    Vector packedMatTransposeTimesVec(const Matrix&, const Vector&);
    void   subtractPackedVecFromVec(Vector& vec, const Vector& packedVec);

    void testAccel(const State&) const;
    void testProjectedVec(const State&, const Vector&) const;
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
    LengthConstraints(const SimbodyMatterSubsystemRep&, int verbose);

    void construct(const Array<RBDistanceConstraint*>&);

    // Returns true if any change was made in the state.
    bool enforcePositionConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;
    bool enforceVelocityConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    // After constraints have been enforced, call these to project out the components
    // of the passed-in q- and u-basis vectors which are normal to the constraint manifold.
    // These project in the Euclidean norm; if you want to project in the error norm
    // you have to scale q and u on the way in and out yourself.
    void projectQVecOntoConfigurationConstraints(const State&, Vector& q);
    void projectUVecOntoMotionConstraints(const State&, Vector& u);

    bool calcConstraintForces(const State&) const;
    void addInCorrectionForces(const State&, SpatialVecList& spatialForces) const;

    void fixVel0(State&, Vector&);


private:
    int    maxIters;
    int    maxMin;

    const SimbodyMatterSubsystemRep& rbTree;
    const int                        verbose;

    std::vector<LengthSet> pvConstraints;   // used for pos, vel
    std::vector<LengthSet> accConstraints;  // used for acc
    NewtonRaphson          posMin, velMin;

    friend class LengthSet;
};

inline const SimbodyMatterSubsystemRep& 
LengthSet::getRBTree()  const {return lConstraints->rbTree;}
inline int                  
LengthSet::getVerbose() const {return lConstraints->verbose;}


#endif // SimTK_SIMBODY_LENGTH_CONSTRAINTS_H_


