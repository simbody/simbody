/**@file
 *
 * Module for bond-length constraints use a shake-like algorithm.
 * Requires integrator to specify time derivatives to constrain.
 * 
 * setup:
 * 1) determine couplings
 * 2) determine which time derivatives to be applied
 * 
 * execution:
 * 1) need all time derivs of theta passed in
 * 2) iteratively apply constraints from lowest to highest time deriv
 * 
 * Deal with loop bond-length constraints.
 */

#include "simbody/internal/common.h"
using namespace simtk;

#include "LengthConstraints.h"
#include "RigidBodyNode.h"
#include "RigidBodyTree.h"

#include "newtonRaphson.h"

#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::setw;

#include <vector>
#include <algorithm>

class LoopWNodes;
static int compareLevel(const LoopWNodes& l1,    //forward declarations
                        const LoopWNodes& l2);

//static bool sameBranch(const RigidBodyNode* tip,
//                       const LoopWNodes& l );

class BadNodeDef {};  //exception

/**
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
    LoopWNodes() : tree(0), rbDistCons(0), flipStations(false), base(0), moleculeNode(0) {}
    LoopWNodes(const RigidBodyTree&, const RBDistanceConstraint&);

    void calcPosInfo(const SBStateRep& s) const { rbDistCons->calcPosInfo(s); }
    void calcVelInfo(const SBStateRep& s) const { rbDistCons->calcVelInfo(s); }
    void calcAccInfo(const SBStateRep& s) const { rbDistCons->calcAccInfo(s); }

    // TODO: State isn't currently in use but will be necessary when the distances
    // are parameters.
    const Real& getDistance(const SBStateRep&) const { return rbDistCons->getDistance(); }

    // Return one of the stations, ordered such that tips(1).level <= tips(2).level.
    const RBStation& tips(int i) const {return rbDistCons->getStation(ix(i));}
    const RigidBodyNode& tipNode(int i) const {return tips(i).getNode();}
 
    const Vec3& tipPos(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return s.getConfigurationCache(*tree).pos_G[ix(i)-1][dc];
    }
    const Vec3& tipVel(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return s.getMotionCache(*tree).vel_G[ix(i)-1][dc];
    }
    const Vec3& tipAcc(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return s.getReactionCache(*tree,true).acc_G[ix(i)-1][dc];
    }
    const Vec3& tipForce(const SBStateRep& s, int i) const {
        assert(1 <= i && i <= 2);
        const int dc = rbDistCons->getDistanceConstraintNum();
        return s.getReactionCache(*tree,true).force_G[ix(i)-1][dc];
    }

    // Use this for both forces and impulses.
    void setTipForce(const SBStateRep& s, int i, const Vec3& f) const {
        assert(1 <= i && i <= 2);

        const int dc = rbDistCons->getDistanceConstraintNum();
        s.updReactionCache(*tree).force_G[ix(i)-1][dc] = f;
    }

private:
    int ix(int i) const { assert(i==1||i==2); return flipStations ? 3-i : i; }
    //const RBStationRuntime& getTipRuntime(int i) const
   //   { return rt->stationRuntimes[ix(i)-1]; }
   // RBStationRuntime& updTipRuntime(int i) const
    //  { return rt->stationRuntimes[ix(i)-1]; }

    const RigidBodyTree*         tree;        // a reference to the tree we're part of
    const RBDistanceConstraint*  rbDistCons;  // reference to the constraint

    // calculated construction-time (topological) info about the constraint
    bool                              flipStations; // make sure station(1).level
                                                    //   <= station(2).level
    std::vector<const RigidBodyNode*> nodes[2];     // the two paths: base..tip1, base..tip2,
                                                    //   incl. tip nodes but not base
    const RigidBodyNode*              base;         // highest-level common ancestor of tips
    const RigidBodyNode*              moleculeNode;

    friend class LengthSet;
    friend ostream& operator<<(ostream& os, const LengthSet& s);
    friend void LengthConstraints::construct(const Array<RBDistanceConstraint*>&);
    friend int compareLevel(const LoopWNodes& l1,
                            const LoopWNodes& l2);
    friend bool sameBranch(const RigidBodyNode* tip,
                           const LoopWNodes& l );
};

LoopWNodes::LoopWNodes(const RigidBodyTree& t, const RBDistanceConstraint& dc)
  : tree(&t), rbDistCons(&dc), flipStations(false), base(0), moleculeNode(0)
{
    const RigidBodyNode* dcNode1 = &dc.getStation(1).getNode();
    const RigidBodyNode* dcNode2 = &dc.getStation(2).getNode();

    if (dcNode1==dcNode2) {
        cout << "LoopWNodes::LoopWNodes: bad topology:\n\t"
             << "loop stations " << dc.getStation(1)
             << " and  "         << dc.getStation(2)
             << " are now in the same node. Deleting loop.\n";
        SIMTK_THROW1(simtk::Exception::LoopConstraintConstructionFailure, "bad topology");
    }

    // Ensure that tips(2) is the atom which is farther from the base.
    flipStations = (dcNode1->getLevel() > dcNode2->getLevel());

    // OK to use tips() at this point.

    // Collect up the node path from tips(2) down to the last node on its
    // side of the loop which is at a higher level than tips(1) (may be none).
    const RigidBodyNode* node1 = &tips(1).getNode();
    const RigidBodyNode* node2 = &tips(2).getNode();
    while ( node2->getLevel() > node1->getLevel() ) {
        nodes[1].push_back(node2);
        node2 = node2->getParent();
    }

    // We're at the same level on both sides of the loop. Run down both
    // simultaneously until we hit the first common ancestor, collecting
    // up the nodes along the two paths (but not the common ancestor).
    while ( node1 != node2 ) {
        if ( node1->isGroundNode() ) {
            cerr << "LoopWNodes::LoopWNodes: could not find base node.\n\t"
                 << "loop between stations " << tips(1) << " and " 
                 << tips(2) << "\n";
            SIMTK_THROW1(simtk::Exception::LoopConstraintConstructionFailure, 
                         "could not find base node");
        }
        nodes[0].push_back(node1);
        nodes[1].push_back(node2);
        node1 = node1->getParent();
        node2 = node2->getParent();
    }

    base = node1;   // that's the common ancestor

    // We want these in base-to-tip order.
    std::reverse(nodes[0].begin(), nodes[0].end());
    std::reverse(nodes[1].begin(), nodes[1].end());

    // find molecule node (level==1 node)
    for (moleculeNode=base; 
         moleculeNode->getLevel()>1 ; 
         moleculeNode=moleculeNode->getParent()) ;

    /* TODO: sherm 060228: allow ground as base node
    if ( moleculeNode->getLevel()<1 ) {
        cerr << "LoopWNodes::LoopWNodes: could not find molecule node.\n\t"
             << "loop between atoms " << tips(1) << " and " 
             << tips(2) << "\n";
        SIMTK_THROW1(simtk::Exception::LoopConstraintConstructionFailure, 
                     "could not find 'molecule' node");
    }
    */
}

typedef std::vector<LoopWNodes> LoopList;

class LengthSet {
    static void construct(const LoopList& loops);
    const LengthConstraints*          lConstraints;
    LoopList                          loops;    
    int                               ndofThisSet;
    std::vector<const RigidBodyNode*> nodeMap; //unique nodes (union of loops->nodes)
public:
    LengthSet() : lConstraints(0), ndofThisSet(0) { }
    LengthSet(const LengthConstraints* lConstraints, const LoopWNodes& loop)
        : lConstraints(lConstraints), ndofThisSet(0)
    {
        addConstraint(loop);
    }

    const RigidBodyTree& getRBTree()  const {return lConstraints->rbTree;}
    int                  getVerbose() const {return lConstraints->verbose;}

    void addConstraint(const LoopWNodes& loop) {
        loops.push_back( LoopWNodes(loop) );
        for (int b=0 ; b<2 ; b++)
            for (int i=0; i<(int)loop.nodes[b].size(); i++)
                if (std::find(nodeMap.begin(),nodeMap.end(),loop.nodes[b][i])
                    ==nodeMap.end())
                {
                    // not found
                    ndofThisSet += loop.nodes[b][i]->getDOF();
                    nodeMap.push_back( loop.nodes[b][i] );
                }
    }

    bool contains(const RigidBodyNode* node) {
        bool found=false;
        for (size_t i=0 ; i<loops.size() ; i++) {
            const std::vector<const RigidBodyNode*>& n0 = loops[i].nodes[0];
            const std::vector<const RigidBodyNode*>& n1 = loops[i].nodes[1];
            if (   (std::find(n0.begin(),n0.end(),node) != n0.end())
                || (std::find(n1.begin(),n1.end(),node) != n1.end()))
            {
                found=true;
                break;
            }
        }
        return found;
    }

    void  setPos(SBStateRep&, const Vector& pos) const;
    void  setVel(SBStateRep&, const Vector& vel) const;
    Vector getPos();
    Vector calcPosB(SBStateRep&, const Vector& pos) const;
    Vector calcVelB(SBStateRep&, const Vector& pos, const Vector& vel) const;
    Vector calcPosZ(const SBStateRep&, const Vector& b) const;
    Matrix calcGrad(const SBStateRep&) const;
    Matrix calcGInverse(const SBStateRep&) const;
    Matrix calcGInverseFD(const SBStateRep&) const;

    void  calcConstraintForces(const SBStateRep&) const; // updates runtime only
    void  addInCorrectionForces(const SBStateRep&, SpatialVecList& spatialForces) const; // spatialForces+=correction

    void   fixVel0(SBStateRep&, Vector&);
    void   fixInternalForce(const SBStateRep&, Vector&);
    Vector multiForce(const Vector&, const Matrix& mat);
    void   subtractVecFromForces(Vector& frc, const Vector& vec);

    void testAccel(const SBStateRep&) const;
    void testInternalForce(const SBStateRep&, const Vector&) const;
    friend ostream& operator<<(ostream& os, const LengthSet& s);
    friend class CalcPosB;
    friend class CalcVelB;
    friend class CalcPosZ;
    friend class CalcVelZ;

    void fdgradf(SBStateRep&, const Vector& pos, Matrix& grad) const;
    void testGrad(SBStateRep&, const Vector& pos, const Matrix& grad) const;
};

ostream& 
operator<<(ostream& os, const LengthSet& s) 
{
    for (size_t i=0 ; i<s.loops.size() ; i++)
        os << setw(4) << s.loops[i].tips(1) << "<->" 
           << setw(4) << s.loops[i].tips(2)  << ": " 
           << s.loops[i].getDistance(SBStateRep()) << "  ";
    return os;
}

class LengthConstraintsPrivates {
public:
    LengthConstraintsPrivates() : posMin(cout), velMin(cout) {}

public:
    std::vector<LengthSet> constraints;     // used for pos, vel
    std::vector<LengthSet> accConstraints;  // used for acc
    NewtonRaphson          posMin, velMin;
};

LengthConstraints::LengthConstraints
    (const RigidBodyTree& rbt, const Real& ctol, int vbose)
  : bandCut(1e-7), maxIters( 20 ), maxMin( 20 ), 
    rbTree(rbt), verbose(vbose), priv(0)
{
    priv = new LengthConstraintsPrivates;
    priv->posMin.maxIters = maxIters;
    priv->posMin.maxMin   = maxMin;
    priv->posMin.tol      = ctol;
    priv->velMin.maxIters = maxIters;
    priv->velMin.maxMin   = maxMin;
    priv->velMin.tol      = ctol*10;
}

LengthConstraints::~LengthConstraints()
{
    delete priv;
}

//
// compare LoopWNodes by base node level value
//
static int
compareLevel(const LoopWNodes& l1,
             const LoopWNodes& l2) 
{ 
    if ( l1.base->getLevel() > l2.base->getLevel() ) 
        return 1;
    else if ( l1.base->getLevel() < l2.base->getLevel() )
        return -1;
    else 
        return 0;
}

// Define this operator so std::sort will work.
static inline bool operator<(const LoopWNodes& l1, const LoopWNodes& l2) {
    return compareLevel(l1,l2) == -1;
}

//1) construct: given list of loops 
//   a) if appropriate issue warning and exit.
//   b) sort by base node of each loop.
//   c) find loops which intersect: combine loops and increment
//    number of length constraints
void
LengthConstraints::construct(const Array<RBDistanceConstraint*>& iloops)
{
    //clean up
    priv->constraints.resize(0);
    priv->accConstraints.resize(0);

    // if ( !useLengthConstraint )  FIX:!!!
    //   //issue error message?
    //   return;

    LoopList loops;
    for (int i=0 ; i < (int)iloops.size() ; i++) {
        try {
            LoopWNodes loop(rbTree, *iloops[i] );
            loops.push_back( loop );
        }
        catch ( BadNodeDef ) {}
    }

    // sort loops by base->level
    std::sort(loops.begin(), loops.end()); // uses "<" operator by default; see above
    //loops.sort(compareLevel);
    LoopList accLoops = loops;  //version for acceleration

    // find intersections -- this version keeps hierarchical loops distinct
    // sherm: the loops are considered coupled if a lower one includes the
    // first body up from the base along either branch of the upper one (if
    // it doesn't include that it can't include any further up the branch either).
    // This makes good sense to me.
    for (int i=0 ; i<(int)loops.size() ; i++) {
        priv->constraints.push_back(LengthSet(this, loops[i]));
        for (int j=i+1 ; j<(int)loops.size() ; j++)
            for (int k=0 ; k<(int)priv->constraints.size() ; k++)
                if (   (loops[j].nodes[0].size()
                        && priv->constraints[k].contains(loops[j].nodes[0][0]))
                    || (loops[j].nodes[1].size()
                        && priv->constraints[k].contains(loops[j].nodes[1][0])))
                {
                    //add length constraint to loop i
                    priv->constraints[i].addConstraint(loops[j]);
                    loops.erase(loops.begin() + j); // STL for &loops[j]
                    j--;
                    break;
                }
    }

    // find intersections - group all loops with tip->trunk relationship
    // TODO sherm: I don't understand why these have to be more coupled than
    // the pos/vel constraints. (This code couples all the loops on the same
    // molecule). I asked Schwieters and he didn't know either but said he
    // would have to think about it. He thought maybe he had cut some corners here
    // and over-coupled the loops but he wasn't sure.
    for (int i=0 ; i<(int)accLoops.size() ; i++) {
        priv->accConstraints.push_back(LengthSet(this, accLoops[i]));
        for (int j=i+1 ; j<(int)accLoops.size() ; j++)
            if ( accLoops[i].moleculeNode == accLoops[j].moleculeNode ) {
                //if (!accLoops[i].moleculeNode->isGroundNode()) { // sherm 060228
                    priv->accConstraints[i].addConstraint(accLoops[j]);
                    accLoops.erase(accLoops.begin() + j); // STL for &accLoops[j]
                    j--;
                //}
            }
        //     for (int b=1 ; b<=2 ; b++) 
        //     if ( sameBranch(accLoops[i].tips(b)->node,accLoops[j]) ||
        //          sameBranch(accLoops[j].tips(b)->node,accLoops[i])   ) {
        //       accConstraints[i]->addConstraint(accLoops[j]);
        //       accLoops.remove(j);
        //       j--;
        //       break;
        //     }
    }

    if (priv->constraints.size()>0 && verbose&InternalDynamics::printLoopInfo) {
        cout << "LengthConstraints::construct: pos/vel length constraints found:\n";
        for (int i=0 ; i<(int)priv->constraints.size() ; i++)
            cout << "\t" << priv->constraints[i] << "\n";
        cout << "LengthConstraints::construct: accel length constraints found:\n";
        for (int i=0 ; i<(int)priv->accConstraints.size() ; i++)
            cout << "\t" << priv->accConstraints[i] << "\n";
    }
}

class CalcPosB {
    SBStateRep&   s;
    const LengthSet* lengthSet;
public:
    CalcPosB(SBStateRep& ss, const LengthSet* lset)
        : s(ss), lengthSet(lset) {}
    Vector operator()(const Vector& pos) const
        { return lengthSet->calcPosB(s,pos); }
};

class CalcPosZ {
    const SBStateRep&   s;
    const LengthSet* lengthSet;
public:
    CalcPosZ(const SBStateRep& ss, const LengthSet* constraint) 
        : s(ss), lengthSet(constraint) {}
    Vector operator()(const Vector& b) const 
        { return lengthSet->calcPosZ(s, b); }
};

class CalcVelB {
    SBStateRep&         s;
    const Vector&    pos;
    const LengthSet* lengthSet;
public:
    CalcVelB(SBStateRep& ss, const Vector& q, const LengthSet* constraint)
        : s(ss), pos(q), lengthSet(constraint) {}
    Vector operator()(const Vector& vel) 
        { return lengthSet->calcVelB(s,pos,vel); }
};

//
// Calculate the position constraint violation (zero when constraint met).
//
Vector
LengthSet::calcPosB(SBStateRep& s, const Vector& pos) const
{
    setPos(s, pos);

    Vector b( loops.size() );
    for (int i=0 ; i<(int)loops.size() ; i++) 
        b(i) = loops[i].getDistance(s) - 
                (loops[i].tipPos(s,1) - loops[i].tipPos(s,2)).norm();
    return b;
}

//
// Calculate the velocity constraint violation (zero when constraint met).
//
Vector
LengthSet::calcVelB(SBStateRep& s, const Vector& pos, const Vector& vel) const 
{
    setPos(s, pos); // TODO (sherm: this is probably redundant)
    setVel(s, vel);

    Vector b( loops.size() );
    for (int i=0 ; i<(int)loops.size() ; i++) {
        // TODO why the minus sign here? (doesn't work right without it) sherm
        b(i) = -dot(unitVec(loops[i].tipPos(s,2) - loops[i].tipPos(s,1)),
                    loops[i].tipVel(s,2) - loops[i].tipVel(s,1));
    }

    return b;
}

//
// Given a vector containing violations of position constraints, calculate
// a state update which should drive those violations to zero (if they
// were linear).
//
// This is a little tricky since the gradient we have is actually the
// gradient of the *velocity* errors with respect to the generalized speeds,
// but we want the gradient of the *position* errors with respect to 
// the generalized coordinates. The theory goes something like this:
//
//       d perr    d perr_dot   d verr    d u
//       ------  = ---------- = ------ * -----
//        d q      d qdot        d u     d qdot
//
// Let P = d perr/dq, G = d verr/du, Q = d qdot/du.
//
// We want to find a change deltaq that will elimate the current error -b:
// P deltaq = b. Instead we solve G * x = b, where x = inv(Q) * deltaq,
// and then solve deltaq = Q * x. Conveniently Q is our friendly invertible
// relation between qdot's and u's: qdot = Q*u.
//
// TODO: I have oversimplified the above since we are really looking for
// a least squares solution to an underdetermined system. With the 
// pseudoinverse G+ available we can write x = G+ b. This is hidden because
// the calcGInverse() routine does return a pseudoinverse, more or less.
//
Vector
LengthSet::calcPosZ(const SBStateRep& s, const Vector& b) const
{
    const Vector x = calcGInverse(s) * b;

    Vector       zu(getRBTree().getTotalDOF(),0.0);
    Vector       zq(getRBTree().getTotalQAlloc(),0.0);

    // map the vector dir back to the appropriate elements of z
    int indx=0; // sherm 060222: I added this
    for (int i=0 ; i<(int)nodeMap.size() ; i++) {
        const int d    = nodeMap[i]->getDOF();
        const int offs = nodeMap[i]->getUIndex();
        zu(offs,d) = x(indx,d);
        indx += d;

        // Make qdot = Q*u.
        nodeMap[i]->calcQDot(s,zu,zq);  // change u's to qdot's
    }
    assert(indx == ndofThisSet);

    return zq;
}

//
// Given a vector containing violations of velocity constraints, calculate
// a state update which would drive those violations to zero (if they
// are linear and well conditioned).
//
// This is simpler than CalcVelP because we have the right gradient here (it's
// the same one in both places).
//
class CalcVelZ {
    const SBStateRep&   s;
    const LengthSet* lengthSet;
    const Matrix     gInverse;
public:
    CalcVelZ(const SBStateRep& ss, const LengthSet* lset)
      : s(ss), lengthSet(lset), gInverse(lengthSet->calcGInverse(s)) {}

    Vector operator()(const Vector& b) {
        const Vector dir = gInverse * b;
        Vector       z(lengthSet->getRBTree().getTotalDOF(),0.0);

        // map the vector dir back to the appropriate elements of z
        int indx=0; // sherm 060222: I added this
        for (int i=0 ; i<(int)lengthSet->nodeMap.size() ; i++) {
            const RigidBodyNode* n = lengthSet->nodeMap[i];
            const int d    = n->getDOF();
            const int offs = n->getUIndex();
            z(offs,d) = dir(indx,d);
            indx += d;
        }
        assert(indx == lengthSet->ndofThisSet);
        return z;
    }
};

// Project out the position constraint errors from the given state. 
bool
LengthConstraints::enforceConfigurationConstraints(SBStateRep& s) const
{
    assert(s.getStage(rbTree) >= ConfiguredStage-1);
    Vector& pos = rbTree.updQ(s);

    bool anyChanges = false;

    try { 
        for (int i=0 ; i<(int)priv->constraints.size() ; i++) {
            anyChanges = true; // TODO: assuming for now
            priv->posMin.calc(pos,
                              CalcPosB(s, &priv->constraints[i]),
                              CalcPosZ(s, &priv->constraints[i]));
        }
    }
    catch ( simtk::Exception::NewtonRaphsonFailure cptn ) {
        cout << "LengthConstraints::enforceConfigurationConstraints: exception: "
             << cptn.getMessage() << '\n';
    } 

    return anyChanges;
}

// Project out the velocity constraint errors from the given state. 
bool
LengthConstraints::enforceMotionConstraints(SBStateRep& s) const
{
    assert(s.getStage(rbTree) >= MovingStage-1);
    const Vector& pos = rbTree.getQ(s);
    Vector&       vel = rbTree.updU(s);

    bool anyChanges = false;

    try { 
        for (int i=0 ; i<(int)priv->constraints.size() ; i++) {
            anyChanges = true; // TODO: assuming for now
            priv->velMin.calc(vel,
                              CalcVelB(s, pos, &priv->constraints[i]),
                              CalcVelZ(s, &priv->constraints[i]));
        }
    }
    catch ( simtk::Exception::NewtonRaphsonFailure cptn ) {
        cout << "LengthConstraints::enforceMotionConstraints: exception: "
             << cptn.getMessage() << '\n';
    } 

    return anyChanges;
}

//
//// on each iteration
//// for each loop
//// -first calc all usual properties
//// -recursively compute phi_ni for each length constraint
////   from tip to base of 
//// -compute gradient
//// -update theta using quasi-Newton-Raphson
//// -compute Cartesian coords
//// -check for convergence
//// -repeat
//


void LengthSet::setPos(SBStateRep& s, const Vector& pos) const
{
    for (int i=0 ; i<(int)nodeMap.size() ; i++)
        nodeMap[i]->setQ(s, pos);

    // TODO: sherm this is the wrong place for the stage update!
    if (s.getStage(getRBTree()) >= ConfiguredStage)
        s.setStage(getRBTree(), TimedStage); // back up if needed

    // sherm TODO: this now computes kinematics for the whole system,
    // but it should only have to update the loop we are interested in.
    // Schwieters had this right before because his equivalent of 'setQ'
    // above also performed the kinematics, while ours just saves the
    // new state variable values and calculates here:
    getRBTree().realizeConfiguration(s);

    // TODO: This is redundant after realizeConfiguration(), but I'm leaving
    // it here because this is actually all that need be recalculated for
    // the loop closure iterations.
    for (int i=0; i<(int)loops.size(); ++i)
        loops[i].calcPosInfo(s);
}

// Must have called LengthSet::setPos() already.
void LengthSet::setVel(SBStateRep& s, const Vector& vel) const
{
    for (int i=0 ; i<(int)nodeMap.size() ; i++)
        nodeMap[i]->setU(s, vel);

    // TODO: sherm this is the wrong place for the stage update!
    if (s.getStage(getRBTree()) >= MovingStage)
        s.setStage(getRBTree(), ConfiguredStage); // back up if needed

    getRBTree().realizeMotion(s);

    // TODO: see comment above in setPos
    for (int i=0; i<(int)loops.size(); ++i)
        loops[i].calcVelInfo(s);
}

//
////A = df / dtheta_i
////  = [ -(q_ni-qn0)x , 1 ] [ phi_na^T Ha , phi_nb^T Hb , ... ]^T
////
//

// Calculate gradient by finite difference for testing the analytic
// version. Presumes that calcEnergy has been called previously with current 
// value of ipos.
void 
LengthSet::fdgradf(SBStateRep& s,
                   const Vector&  pos,
                   Matrix&        grad) const 
{
    // Gradf gradf(tree);
    // gradf(x,grad); return;
    const double eps = 1e-8;
    const CalcPosB calcB(s, this);
    const Vector b = calcB(pos);
    int grad_indx=0;
    for (int i=0 ; i<(int)nodeMap.size() ; i++) {
        int pos_indx=nodeMap[i]->getQIndex();
        for (int j=0 ; j<nodeMap[i]->getNQ(s) ; j++,pos_indx++,grad_indx++) {
            Vector posp = pos;
            posp(pos_indx) += eps;
            const Vector bp = calcB(posp);
            posp(pos_indx) -= 2.*eps;
            const Vector bpm = calcB(posp);
            for (int k=0 ; k<b.size() ; k++)
                grad(grad_indx,k) = -(bp(k)-bpm(k)) / (2.*eps);
        }
    }
}

void 
LengthSet::testGrad(SBStateRep& s, const Vector& pos, const Matrix& grad) const
{
    double tol = 1e-4;

    Matrix fdgrad(ndofThisSet,loops.size());
    fdgradf(s,pos,fdgrad);

    for (int i=0 ; i<grad.nrow() ; i++)
        for (int j=0 ; j<grad.ncol() ; j++)
            if (fabs(grad(i,j)-fdgrad(i,j)) > fabs(tol))
                cout << "testGrad: error in gradient: " 
                     << setw(2) << i << ' '
                     << setw(2) << j << ": "
                     << grad(i,j) << ' ' << fdgrad(i,j) << '\n';
    cout.flush();
}

//
// unitVec(p+ - p-) * d (p+ - p-) / d (theta_i)
// d g / d theta for all hingenodes in nodemap
//
// sherm: this appears to calculate the transpose of G
//
// TODO: won't work right for balls right now
// sherm 060222: OK, this routine doesn't really do what it says
// when the constraint includes a ball or free joint. It
// does not use the actual q's, which are quaternions. Instead
// it is something like d(v+ - v-)/du. If the u's are the 
// angular velocity across the joint then this is
// d(p+ - p-)/d qbar where qbar is a 1-2-3 Euler sequence
// which is 0,0,0 at the current orientation. ("instant
// coordinates").
//
// sherm: 060303 This routine uses the joint transition matrices H,
// which can be thought of as Jacobians d V_PB_G / d uB, that is
// partial derivative of the cross-joint relative *spatial* velocity
// with respect to that joint's generalized speeds uB. This allows
// analytic computation of d verr / d u where verr is the set of
// distance constraint velocity errors for the current set of 
// coupled loops, and u are the generalized speeds for all joints
// which contribute to any of those loops. Some times this is the
// gradient we want, but we also want to use this routine to
// calculate d perr / d q, which we can't get directly. But it
// is easily finagled into the right gradient; see CalcPosZ above.
// The trickiest part is that we use nice physical variables for
// generalized speeds, like angular velocities for ball joints,
// but we use awkward mathematical constructs like quaternions and
// Euler angles for q's.
//
Matrix
LengthSet::calcGrad(const SBStateRep& s) const
{
    Matrix grad(ndofThisSet,loops.size(),0.0);
    const Mat33 one(1);  //FIX: should be done once

    for (int i=0 ; i<(int)loops.size() ; i++) {
        const LoopWNodes& l = loops[i];
        std::vector<SpatialMat> phiT[2];
        for (int b=0 ; b<2 ; b++) {
            phiT[b].resize( l.nodes[b].size() );
            if ( l.nodes[b].size() ) {
                phiT[b][phiT[b].size()-1] = 1;  // identity
                for (int j=l.nodes[b].size()-2 ; j>=0 ; j-- ) {
                    const RigidBodyNode* n = l.nodes[b][j+1];
                    phiT[b][j] = phiT[b][j+1] * ~n->getPhi(s);
                }
            }
        }

        // compute gradient
        Vec3 uBond = unitVec(l.tipPos(s,2) - l.tipPos(s,1));
        Row<2,Mat33> J[2];
        for (int b=1 ; b<=2 ; b++)
            // TODO: get rid of this b-1; make tips 0-based
            J[b-1] = Row<2,Mat33>(-crossMat(l.tipPos(s,b) -
                                            l.tips(b).getNode().getX_GB(s).T()),   one);
        int g_indx=0;
        for (int j=0 ; j<(int)nodeMap.size() ; j++) {
            Real elem=0.0;

            // We just want to get the index at which nodeMap[j] is found in the
            // std::vector (or -1 if not found) but that's not so easy!
            const std::vector<const RigidBodyNode*>& n0 = l.nodes[0];
            const std::vector<const RigidBodyNode*>& n1 = l.nodes[1];
            std::vector<const RigidBodyNode*>::const_iterator found0 =
                std::find(n0.begin(),n0.end(),nodeMap[j]);
            std::vector<const RigidBodyNode*>::const_iterator found1 =
                std::find(n1.begin(),n1.end(),nodeMap[j]);

            const int l1_indx = (found0==n0.end() ? -1 : found0-n0.begin());
            const int l2_indx = (found1==n1.end() ? -1 : found1-n1.begin());

            for (int k=0 ; k < nodeMap[j]->getDOF() ; k++) {
                const SpatialVec& HtCol = ~nodeMap[j]->getHRow(s, k);
                if ( l1_indx >= 0 ) { 
                    elem = -dot(uBond , Vec3(J[0] * phiT[0][l1_indx]*HtCol));
                } else if ( l2_indx >= 0 ) { 
                    elem =  dot(uBond , Vec3(J[1] * phiT[1][l2_indx]*HtCol));
                }
                grad(g_indx++,i) = elem;
            }
        }
        // added by sherm 060222:
        assert(g_indx == ndofThisSet); // ??
    }
    return grad;
} 

//
// Calculate generalized inverse which minimizes changes in soln vector.
// TODO (sherm) This is trying to create a pseudoinverse
// using normal equations which is numerically bad and can't deal with
// redundant constraints. Should use an SVD or (faster) QTZ factorization 
// instead.
//
// sherm 060314:
//                           n
//                  --------------------
//                 |                    |
//     A (mXn) = m |                    |   rank(A) <= m.
//                 |                    |
//                  --------------------
//
// The nXm pseudo-inverse A+ of a matrix A is A+ = V S+ ~U, where A=U*S*~V
// is the singular value decomposition (SVD) of A, and S+ is the inverse
// of the diagonal matrix S with some zeroes thrown it at the end after
// we run out of rank(A). For an underdetermined system with full row
// rank, that is, assuming m < n and rank(A)=m, you can do a poor man's
// calculation of this as A+ = ~A*inv(A*~A). In the overdetermined
// case we have m > n and rank(A)=n, in which case A+ = inv(~A*A)*~A.
//
// We are given gradient G (nuXnc). We're assuming nu > nc and rank(G)=nc, i.e.,
// no redundant constraints (not a good assumption in general!). We would like
// to get a least squares solution to ~G x = b where x has dimension nu and
// b has dimension nc. So if we set A=~G we have the underdetermined system
// depicted above and want to return A+ = ~A*inv(A*~A) = G*inv(~G*G). Ergo ...
Matrix
LengthSet::calcGInverse(const SBStateRep& s) const
{
    const Matrix A = ~calcGrad(s); // now A is dg/dtheta
    const int m = A.nrow(); // see picture above
    const int n = A.ncol();

    // Now calculate the pseudo inverse of A.
    Matrix pinvA(n,m,0.);
    if ( A.normSqr() > 1e-10 ) {
        if (m < n)
            pinvA = ~A * (A * ~A).invert(); // normal underdetermined case
        else if (m > n)
            pinvA = (A * ~A).invert() * ~A; // wow, that's a lot of constraints!
        else 
            pinvA = A.invert();             // not likely
    }
    return pinvA;
}

Matrix
LengthSet::calcGInverseFD(const SBStateRep& s) const
{
    Matrix grad(ndofThisSet,loops.size()); // <-- transpose of the actual dg/dtheta
    SBStateRep sTmp = s;
    Vector pos = getRBTree().getQ(s);
    fdgradf(sTmp,pos,grad);

    const Matrix A = ~grad; // now A is dg/dtheta
    const int m = A.nrow(); // see picture above
    const int n = A.ncol();

    // Now calculate the pseudo inverse of A.
    Matrix pinvA(n,m,0.);
    if ( A.normSqr() > 1e-10 ) {
        if (m < n)
            pinvA = ~A * (A * ~A).invert(); // normal underdetermined case
        else if (m > n)
            pinvA = (A * ~A).invert() * ~A; // wow, that's a lot of constraints!
        else 
            pinvA = A.invert();             // not likely
    }
    return pinvA;
}

//acceleration:
//  0) after initial acceleration calculation:
//  1) calculate Y (block diagonal matrix) for all nodes (base to tip)
//  2) for each LengthSet, calculate Lagrange multiplier(s) and add in
//     resulting force- update accelerations
//     Do this step from tip to base.
//

bool LengthConstraints::calcConstraintForces(const SBStateRep& s) const {
    if ( priv->accConstraints.size() == 0 )
        return false;

    rbTree.calcY(s);

    for (int i=priv->accConstraints.size()-1 ; i>=0 ; i--)
        priv->accConstraints[i].calcConstraintForces(s);

    return true;
}

void LengthConstraints::addInCorrectionForces(const SBStateRep& s, SpatialVecList& spatialForces) const {
    for (int i=priv->accConstraints.size()-1 ; i>=0 ; i--)
        priv->accConstraints[i].addInCorrectionForces(s, spatialForces);
}

void 
LengthConstraints::fixGradient(const SBStateRep& s, Vector& forceInternal)
{
    if ( priv->constraints.size() == 0 )
        return;

    for (int i=priv->constraints.size()-1 ; i>=0 ; i--)
        priv->constraints[i].fixInternalForce(s, forceInternal);

    for (int i=priv->constraints.size()-1 ; i>=0 ; i--) 
        priv->constraints[i].testInternalForce(s, forceInternal);
}

//
// Calculate the acceleration of atom assuming that the spatial acceleration
// of the body (node) it's on is available.
//
//static Vec3
//getAccel(const IVMAtom* a)
//{
//    const RigidBodyNode* n = a->node;
//    Vec3 ret( SubVec6(n->getSpatialAcc(),3,3).vector() );
//    ret += cross( SubVec6(n->getSpatialAcc(),0,3).vector() , a->pos - n->getAtom(0)->pos );
//    ret += cross( SubVec6(n->getSpatialVel(),0,3).vector() , a->vel - n->getAtom(0)->vel );
//    return ret;
//}

//             T     -1 T   
// Compute   v1 *(J M  J )_mn * v2   where the indices mn are given by 
// the nodes associated with stations s1 & s2.
//
static double
computeA(const SBStateRep& s,
         const Vec3&    v1,
         const LoopWNodes& loop1, int s1,
         const LoopWNodes& loop2, int s2,
         const Vec3&    v2)
{
    const RigidBodyNode* n1 = &loop1.tips(s1).getNode();
    const RigidBodyNode* n2 = &loop2.tips(s2).getNode();

    const Mat33 one(1);

    SpatialRow t1 = ~v1 * Row<2,Mat33>(crossMat(n1->getX_GB(s).T() - loop1.tipPos(s,s1)), one);
    SpatialVec t2 = Vec<2,Mat33>(crossMat(loop2.tipPos(s,s2) - n2->getX_GB(s).T()), one) * v2;

    while ( n1->getLevel() > n2->getLevel() ) {
        t1 = t1 * ~n1->getPsi(s);
        n1 = n1->getParent();
    }

    while ( n2->getLevel() > n1->getLevel() ) {
        t2 = n2->getPsi(s) * t2;
        n2 = n2->getParent();
    }

    while ( n1 != n2 ) {
        if (n1->isGroundNode() || n2->isGroundNode()  ) {
            t1 = 0.;  //not in same branch (or same tree -- sherm)
            t2 = 0.;
            cout << "computeA: cycles wasted calculating missed branch: "
                    << loop1.tips(s1) << " <-> " << loop2.tips(s2) << '\n';
            return 0.;
        }
        t1 = t1 * ~n1->getPsi(s);
        t2 = n2->getPsi(s) * t2;
        n1 = n1->getParent();
        n2 = n2->getParent();
    }

    // here n1==n2

    double ret = t1 * n1->getY(s) * t2;
    return ret;
}

//
// To be called for LengthSets consecutively from tip to base.
// This will calculate a force for every station in every loop
// contained in this LengthSet, and store that force in the runtime
// block associated with that loop. It is up to the caller to do
// something with these forces.
//
// See Section 2.6 on p. 294 of Schwieters & Clore, 
// J. Magnetic Resonance 152:288-302. Equation reference below
// are to that paper. (sherm)
//
void
LengthSet::calcConstraintForces(const SBStateRep& s) const
{
    // This is the acceleration error for each loop constraint in this
    // LengthSet. We get a single scalar error per loop, since each
    // contains one distance constraint.
    // See Eq. [53] and the last term of Eq. [66].
    Vector rhs(loops.size(),0.);
    for (int i=0 ; i<(int)loops.size() ; i++) {
        rhs(i) = (loops[i].tipVel(s,2) - loops[i].tipVel(s,1)).normSqr()
                   + dot(loops[i].tipAcc(s,2) - loops[i].tipAcc(s,1) , 
                         loops[i].tipPos(s,2) - loops[i].tipPos(s,1));
    }

    // Here A = Q*(J inv(M) J')*Q' where J is the kinematic Jacobian for
    // the constrained points and Q is the constraint Jacobian. See first 
    // term of Eq. [66].
    Matrix A(loops.size(),loops.size(),0.);
    for (int i=0 ; i<(int)loops.size() ; i++) {
        const Vec3 v1 = loops[i].tipPos(s,2) - loops[i].tipPos(s,1);
        for (int bi=1 ; bi<=2 ; bi++)
            for (int bj=1 ; bj<=2 ; bj++) {
                double maxElem = 0.;
                for (int j=i ; j<(int)loops.size() ; j++) {
                    const Vec3 v2 = loops[j].tipPos(s,2) - loops[j].tipPos(s,1);
                    double contrib = computeA(s, v1, loops[i], bi,
                                                     loops[j], bj, v2);
                    A(i,j) += contrib * (bi==bj ? 1 : -1);

                    if ( fabs(contrib) > maxElem ) maxElem = fabs(contrib);
                    if ( maxElem>0. && fabs(contrib)/maxElem < lConstraints->bandCut ) {
                        if ( getVerbose()&InternalDynamics::printLoopDebug )
                            cout << "LengthSet::calcConstraintForces: setting A("
                                 << i << "," << j+1 << ".." << loops.size()-1 << ")["
                                 << bi << "," << bj << "] = 0\n";
                        break;  //don't compute smaller elements
                    }
                }
            }
    }

    // (sherm) Ouch -- this part is very crude. If this is known to be well 
    // conditioned it should be factored with a symmetric method
    // like Cholesky, otherwise use an SVD (or better, complete orthogonal
    // factorization QTZ) to get appropriate least squares solution.

    for (int i=0 ; i<(int)loops.size() ; i++)   //fill lower triangle
        for (int j=0 ; j<i ; j++)
            A(i,j) = A(j,i);

    //cout << "Solve A lambda=rhs; A=" << A; 
    //cout << "  rhs = " << rhs << endl;

    //FIX: using inverse is inefficient
    const Vector lambda = A.invert() * rhs;

    //cout << "  lambda = " << lambda << endl;

    // add forces due to these constraints
    for (int i=0 ; i<(int)loops.size() ; i++) {
        const Vec3 frc = lambda(i) * (loops[i].tipPos(s,2) - loops[i].tipPos(s,1));
        loops[i].setTipForce(s, 2, -frc);
        loops[i].setTipForce(s, 1,  frc);
    }
}

void LengthSet::addInCorrectionForces(const SBStateRep& s, SpatialVecList& spatialForces) const {
    for (int i=0; i<(int)loops.size(); ++i) {
        for (int t=1; t<=2; ++t) {
            const RigidBodyNode& node = loops[i].tipNode(t);
            const Vec3 force = loops[i].tipForce(s,t);
            const Vec3 moment = cross(loops[i].tipPos(s,t) - node.getX_GB(s).T(), force);
            spatialForces[node.getNodeNum()] += SpatialVec(moment, force);
        }
    }
}

void LengthSet::testAccel(const SBStateRep& s) const
{
    double testTol=1e-8;
    for (int i=0 ; i<(int)loops.size() ; i++) {
        double test=   dot(loops[i].tipAcc(s,2) - loops[i].tipAcc(s,1),
                           loops[i].tipPos(s,2) - loops[i].tipPos(s,1))
                     + (loops[i].tipVel(s,2)-loops[i].tipVel(s,1)).normSqr();
        if ( fabs(test) > testTol )
            cout << "LengthSet::testAccel: constraint condition between atoms "
                 << loops[i].tips(1) << " and " << loops[i].tips(2) << " violated.\n"
                 << "\tnorm of violation: " << fabs(test) << '\n';
    }
    cout.flush();
}
   
//
// Fix internal force such that it obeys the constraint conditions.
//
void
LengthSet::fixInternalForce(const SBStateRep& s, Vector& forceInternal)
{
    Matrix  grad = calcGrad(s);
    if ( getVerbose() & InternalDynamics::printLoopDebug ) {
        SBStateRep sTmp = s;
        Vector pos = getRBTree().getQ(s);
        testGrad(sTmp,pos,grad);
    }
    const Vector rhs  = multiForce(forceInternal,grad); 

    const Matrix A = ~grad * grad;

    //FIX: using inverse is inefficient
    const Vector lambda = A.invert() * rhs;

    // subtract forces due to these constraints
    subtractVecFromForces(forceInternal, grad * lambda); 
}

// This just computes ~M*f but first plucks out the relevant entries
// in f to squash it down to the same size as M.
Vector 
LengthSet::multiForce(const Vector& forceInternal, const Matrix& mat)
{
    // sherm 060222: 
    assert(forceInternal.size() == getRBTree().getTotalDOF());
    assert(mat.nrow() == ndofThisSet && mat.ncol() == ndofThisSet);

    Vector vec(ndofThisSet);    //build vector same size as mat
    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        vec(indx,d) = forceInternal(offs,d);
        indx += d;
    }
    // sherm 060222: 
    assert(indx == ndofThisSet);

    const Vector ret = ~mat * vec; 
    return ret;
}

void
LengthSet::subtractVecFromForces(Vector& forceInternal,
                                 const Vector& vec)
{
    // sherm 060222:
    assert(forceInternal.size() == getRBTree().getTotalDOF());
    assert(vec.size() == ndofThisSet);

    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        forceInternal(offs,d) -= vec(indx,d);
        indx += d;
    }

    assert(indx == ndofThisSet);
}

void
LengthSet::testInternalForce(const SBStateRep& s, 
                             const Vector& forceInternal) const
{
    // sherm 060222:
    assert(forceInternal.size() == getRBTree().getTotalDOF());

    Vector vec(ndofThisSet);
    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        vec(indx,d) = forceInternal(offs,d);
        indx += d;
    }
    assert(indx == ndofThisSet);

    const Matrix grad = calcGrad(s);
    const Vector test = ~grad * vec;

    double testTol=1e-8;
    for (int i=0 ; i<(int)loops.size() ; i++) {
        if ( test(i) > testTol )
            cout << "LengthSet::Gradient: constraint condition between atoms "
                 << loops[i].tips(1) << " and " << loops[i].tips(2) << " violated.\n"
                 << "\tnorm of violation: " << test(i) << '\n';
    }
    cout.flush();
}

void
LengthConstraints::fixVel0(SBStateRep& s, Vector& iVel)
{
    // use molecule grouping
    for (int i=0 ; i<(int)priv->accConstraints.size() ; i++)
        priv->accConstraints[i].fixVel0(s, iVel);
}

//
// Correct internal velocities such that the constraint conditions are
// met under the condition that the disturbance to station velocites is
// small as possible.
//
void
LengthSet::fixVel0(SBStateRep& s, Vector& iVel)
{
    assert(iVel.size() == getRBTree().getTotalDOF());

    // store internal velocities
    Vector iVel0 = iVel;

    // verr stores the current velocity errors, which we're assuming are valid.
    Vector verr(loops.size());
    for (int k=0 ; k<(int)loops.size() ; k++) 
        verr[k] = dot(loops[k].tipPos(s,2) - loops[k].tipPos(s,1) ,
                      loops[k].tipVel(s,2) - loops[k].tipVel(s,1));

    Matrix mat(loops.size(),loops.size());
    std::vector<Vector> deltaIVel(loops.size());
    for (int m=0 ; m<(int)loops.size() ; m++) {
        deltaIVel[m].resize(iVel.size());

        // Set all velocities to zero. TODO: this should just be an "ignore velocity"
        // option to realizeMotion(); it shouldn't actually require putting zeroes everywhere.
        iVel = 0.;
        getRBTree().setU(s, iVel );
        getRBTree().realizeMotion(s);

        // sherm: I think the following is a unit "probe" velocity, projected
        // along the separation vector. 
        // That would explain the fact that there are no velocities here!
        const Vec3 probeImpulse = loops[m].tipPos(s,2)-loops[m].tipPos(s,1);
        loops[m].setTipForce(s, 2,  probeImpulse);
        loops[m].setTipForce(s, 1, -probeImpulse);

        // Convert the probe impulses at the stations to spatial impulses at
        // the body origin.
        SpatialVecList spatialImpulse(getRBTree().getNBodies());
        spatialImpulse.setToZero();
        for (int t=1; t<=2; ++t) {
            const RigidBodyNode& node = loops[m].tipNode(t);
            const Vec3 force = loops[m].tipForce(s,t);
            const Vec3 moment = cross(loops[m].tipPos(s,t) - node.getX_GB(s).T(), force);
            spatialImpulse[node.getNodeNum()] += SpatialVec(moment, force);
        }

        // Apply probe impulse as though it were a force; the resulting "acceleration"
        // is actually the deltaV produced by this impulse, that is, a deltaV which
        // results in a change in velocity along the line between the stations.

        // calc deltaVa from k-th constraint condition
        //FIX: make sure that nodeMap is ordered by level
        // get all nodes in given molecule, ordered by level
        getRBTree().calcZ(s, spatialImpulse);
        getRBTree().calcTreeAccel(s);
        deltaIVel[m] = getRBTree().getUDot(s);

        // Set the new velocities.
        getRBTree().setU(s, deltaIVel[m] );
        getRBTree().realizeMotion(s);

        // Calculating partial(velocityError[n])/partial(deltav[m]). Any velocity
        // we see here is due to the deltav, since we started out at zero.
        for (int n=0 ; n<(int)loops.size() ; n++)
            mat(n,m) = dot(loops[n].tipPos(s,2) - loops[n].tipPos(s,1),
                           loops[n].tipVel(s,2) - loops[n].tipVel(s,1));

        //store results of m-th constraint on deltaVa-n
    }

    // Calculate the fraction of each deltaV by which we should change V to simultaneously
    // drive all the velocity errors to zero.
    // TODO: deal with a badly conditioned Jacobian to produce a least squares lambda. 
    const Vector lambda = mat.invert() * verr;

    iVel = iVel0;
    for (int m=0 ; m<(int)loops.size() ; m++)
        iVel -= lambda[m] * deltaIVel[m];
    getRBTree().setU(s, iVel);
    getRBTree().realizeMotion(s);
}

    // RBDistanceConstraint methods

void RBDistanceConstraint::calcStationPosInfo(const SBStateRep& s, int i) const {
    updStation_G(s,i) = getNode(i).getX_GB(s).R() * getPoint(i);
    updPos_G(s,i)     = getNode(i).getX_GB(s).T() + getStation_G(s,i);
}

void RBDistanceConstraint::calcStationVelInfo(const SBStateRep& s, int i) const {
    const Vec3& w_G = getNode(i).getSpatialAngVel(s);
    const Vec3& v_G = getNode(i).getSpatialLinVel(s);
    updStationVel_G(s,i) = cross(w_G, getStation_G(s,i));
    updVel_G(s,i)        = v_G + getStationVel_G(s,i);
}

void RBDistanceConstraint::calcStationAccInfo(const SBStateRep& s, int i) const {
    const Vec3& w_G  = getNode(i).getSpatialAngVel(s);
    const Vec3& v_G  = getNode(i).getSpatialLinVel(s);
    const Vec3& aa_G = getNode(i).getSpatialAngAcc(s);
    const Vec3& a_G  = getNode(i).getSpatialLinAcc(s);
    updAcc_G(s,i) = a_G + cross(aa_G, getStation_G(s,i))
                        + cross(w_G,  getStationVel_G(s,i)); // i.e., w X (wXr)
}

std::ostream& operator<<(std::ostream& o, const RBStation& s) {
    o << "station " << s.getPoint() << " on node " << s.getNode().getNodeNum();
    return o;
}

void RBDistanceConstraint::calcPosInfo(const SBStateRep& s) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationPosInfo(s, i);

    updFromTip1ToTip2_G(s)  = getPos_G(s,2) - getPos_G(s,1);
    const double separation = getFromTip1ToTip2_G(s).norm();
    updUnitDirection_G(s)   = getFromTip1ToTip2_G(s) / separation;
    updPosErr(s) = distance - separation;
}

void RBDistanceConstraint::calcVelInfo(const SBStateRep& s) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationVelInfo(s, i);

    updRelVel_G(s) = getVel_G(s,2) - getVel_G(s,1);
    updVelErr(s)   = ~getUnitDirection_G(s) * getRelVel_G(s);
}

void RBDistanceConstraint::calcAccInfo(const SBStateRep& s) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationAccInfo(s, i);

//XXX this doesn't look right
    const Vec3 relAcc_G = getAcc_G(s,2) - getAcc_G(s,1);
    updAccErr(s) = getRelVel_G(s).normSqr() + (~relAcc_G * getFromTip1ToTip2_G(s));
}

