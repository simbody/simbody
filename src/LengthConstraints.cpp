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
 * Contributors: Derived from NIH IVM code written by Charles Schwieters      *
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
using namespace SimTK;

#include "LengthConstraints.h"
#include "RigidBodyNode.h"
#include "SimbodyMatterSubsystemRep.h"

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


LoopWNodes::LoopWNodes(const SimbodyMatterSubsystemRep& t, const RBDistanceConstraint& dc)
  : tree(&t), rbDistCons(&dc), flipStations(false), outmostCommonBody(0)
{
    const RigidBodyNode* dcNode1 = &dc.getStation(1).getNode();
    const RigidBodyNode* dcNode2 = &dc.getStation(2).getNode();

    if (dcNode1==dcNode2) {
        cout << "LoopWNodes::LoopWNodes: bad topology:\n\t"
             << "loop stations " << dc.getStation(1)
             << " and  "         << dc.getStation(2)
             << " are now in the same node. Deleting loop.\n";
        SimTK_THROW1(SimTK::Exception::LoopConstraintConstructionFailure, "bad topology");
    }

    // Ensure that tips(2) is on the body which is farther from Ground.
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
            SimTK_THROW1(SimTK::Exception::LoopConstraintConstructionFailure, 
                         "could not find base node");
        }
        nodes[0].push_back(node1);
        nodes[1].push_back(node2);
        node1 = node1->getParent();
        node2 = node2->getParent();
    }

    outmostCommonBody = node1;   // that's the common ancestor; might be Ground

    // We want these in base-to-tip order.
    std::reverse(nodes[0].begin(), nodes[0].end());
    std::reverse(nodes[1].begin(), nodes[1].end());

    // Make the list of ancestors, starting with the outmostCommonBody unless it
    // is ground. Put in base-to-tip order.
    const RigidBodyNode* next = outmostCommonBody;
    while (next->getLevel() != 0) {
        ancestors.push_back(next);
        next = next->getParent();
    }
    std::reverse(ancestors.begin(), ancestors.end());
}


ostream& 
operator<<(ostream& o, const LengthSet& s) 
{
    o << "LengthSet --------------->\n";
    for (int i=0 ; i<(int)s.loops.size() ; i++) {
        o << "Loop " << i << ":\n" << s.loops[i];
    }
    o << "\nNodemap: ";
    for (int i=0; i<(int)s.nodeMap.size(); ++i)
        o << " " << s.nodeMap[i]->getNodeNum()
            << "[" << s.nodeMap[i]->getLevel() << "]";
    o << "\n<-------- end of LengthSet\n";
    return o;
}


ostream& 
operator<<(ostream& o, const LoopWNodes& w) 
{
    o << "tip1=" << w.tips(1) << " tip2=" << w.tips(2) 
      << " distance=" << w.getDistance() << endl;
    o << "nodes[0]:";
    for (int i=0; i<(int)w.nodes[0].size(); ++i)
        o << " " << w.nodes[0][i]->getNodeNum() 
                 << "[" << w.nodes[0][i]->getLevel() << "]";
    o << "\nnodes[1]:";
    for (int i=0; i<(int)w.nodes[1].size(); ++i)
        o << " " << w.nodes[1][i]->getNodeNum() 
                 << "[" << w.nodes[1][i]->getLevel() << "]";
    o << "\nOutmost common body: " << w.getOutmostCommonBody()->getNodeNum()
        << "[" << w.getOutmostCommonBody()->getLevel() << "]";
    o << "\nAncestors: ";
    for (int i=0; i<(int)w.ancestors.size(); ++i)
        o << " " << w.ancestors[i]->getNodeNum()
            << "[" << w.ancestors[i]->getLevel() << "]";
    o << endl;
    return o;
}

LengthConstraints::LengthConstraints
    (const SimbodyMatterSubsystemRep& rbt, int vbose)
  : maxIters( 20 ), maxMin( 20 ), 
    rbTree(rbt), verbose(vbose), posMin("posMin", cout), velMin("velMin", cout)
{
    posMin.maxIters = maxIters;
    posMin.maxMin   = maxMin;
    velMin.maxIters = maxIters;
    velMin.maxMin   = maxMin;
}

//
// compare LoopWNodes by outmost common node level value
//
static int
compareLevel(const LoopWNodes& l1,
             const LoopWNodes& l2) 
{ 
    if ( l1.getOutmostCommonBody()->getLevel() > l2.getOutmostCommonBody()->getLevel() ) 
        return 1;
    else if ( l1.getOutmostCommonBody()->getLevel() < l2.getOutmostCommonBody()->getLevel() )
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
//
// Sherm's constraint coupling hypothesis (060616):
//
// For a constraint i:
//   Kinematic participants K[i] 
//       = the set of bodies whose mobilities can affect verr[i]
//   Dynamic participants D[i]
//       = the set of bodies whose effective inertias can change 
//         as a result of enforcement of constraint i
// TODO: A thought on dynamic coupling -- is this the same as computing
//       effective joint forces? Then only mobilities that can feel a
//       test constraint force can be considered dynamic participants of
//       that constraint.
// 
// Constraints i & j are directly kinematically coupled if 
//   K[i] intersect K[j] <> { }
// Constraints i & j are directly dynamically coupled if
//   (a) they are directly kinematically coupled, or
//   (b) D[i] intersect D[j] <> { }
// TODO: I *think* dynamic coupling is one-way, from outboard to inboard.
//       That means it defines an evaluation *order*, where modified 
//       inertias from outboard loops become the articulated body
//       inertias to use in the inboard loops, rather than requiring simultaneous
//       solutions for the multipliers.
//       Kinematic coupling, on the other hand, is symmetric and implies
//       simultaneous solution.
//
// Compute transitive closure:
//   Constraints i & j are kinematically coupled if
//   (a) they are directly kinematically coupled, or
//   (b) constraint i is kinematically coupled to any constraint
//          k to which j is kinematically coupled
// That is, build kinematic clusters with completely disjoint constraints.
//
// For now, do the same thing with acceleration constraints but this
// is too strict. TODO
// 
//   
void
LengthConstraints::construct(const std::vector<RBDistanceConstraint*>& iloops)
{
    //clean up
    pvConstraints.resize(0);
    accConstraints.resize(0);

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

    //cout << "RAW LOOPS:\n";
    //for (int i=0; i < (int)loops.size(); ++i)
    //    cout << "Loop " << i << ":\n  " << loops[i] << endl;

    // sort loops by outmostCommonBody->level
    std::sort(loops.begin(), loops.end()); // uses "<" operator by default; see above
    LoopList accLoops = loops;  //version for acceleration

    // Sherm 060616: transitive closure calculation. Remember that 
    // the constraints are sorted by the level of their outmost common body, 
    // starting at ground, so that more-outboard constraints cannot cause
    // us to need to revisit an earlier cluster. However, constraints at
    // the same level can bring together earlier ones at that level.
    // TODO: currently just restarting completely when we add a
    // constraint to a cluster.

    for (int i=0 ; i<(int)loops.size() ; i++) {
        pvConstraints.push_back(LengthSet(this));
        pvConstraints.back().addKinematicConstraint(loops[i]);
        bool addedAConstraint;
        do {
            addedAConstraint = false;
            for (int j=i+1 ; j<(int)loops.size() ; j++)
                if (   (loops[j].nodes[0].size()
                        && pvConstraints[i].contains(loops[j].nodes[0][0]))
                    || (loops[j].nodes[1].size()
                        && pvConstraints[i].contains(loops[j].nodes[1][0])))
                {
                    //add constraint equation j to cluster i
                    pvConstraints[i].addKinematicConstraint(loops[j]);
                    loops.erase(loops.begin() + j); // STL for &loops[j]
                    addedAConstraint = true;
                    break;
                }
        } while (addedAConstraint);
    }

    for (int i=0 ; i<(int)accLoops.size() ; i++) {
        accConstraints.push_back(LengthSet(this));
        accConstraints.back().addDynamicConstraint(accLoops[i]);
        bool addedAConstraint;
        do {
            addedAConstraint = false;
            for (int j=i+1 ; j<(int)accLoops.size() ; j++)
                if (   (accLoops[j].nodes[0].size()
                        && accConstraints[i].contains(accLoops[j].nodes[0][0]))
                    || (accLoops[j].nodes[1].size()
                        && accConstraints[i].contains(accLoops[j].nodes[1][0]))
                    || (accLoops[j].ancestors.size()
                        && accConstraints[i].contains(accLoops[j].ancestors[0])))
                {
                    //add constraint equation j to cluster i
                    accConstraints[i].addDynamicConstraint(accLoops[j]);
                    accLoops.erase(accLoops.begin() + j); // STL for &accLoops[j]
                    addedAConstraint = true;
                    break;
                }
        } while (addedAConstraint);
    }

    if (false && pvConstraints.size()>0) {
        cout << "LengthConstraints::construct: pos/vel length sets found:\n";
        for (int i=0 ; i<(int)pvConstraints.size() ; i++)
            cout << pvConstraints[i] << "\n";
        cout << "LengthConstraints::construct: accel length sets found:\n";
        for (int i=0 ; i<(int)accConstraints.size() ; i++)
            cout << accConstraints[i] << "\n";
    }
}

void LengthSet::addKinematicConstraint(const LoopWNodes& loop) {
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

void LengthSet::addDynamicConstraint(const LoopWNodes& loop) {
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

    for (int i=0; i<(int)loop.ancestors.size(); i++)
        if (std::find(nodeMap.begin(),nodeMap.end(),loop.ancestors[i])
            ==nodeMap.end())
        {
            // not found
            ndofThisSet += loop.ancestors[i]->getDOF();
            nodeMap.push_back( loop.ancestors[i] );
        }
}

bool LengthSet::contains(const RigidBodyNode* node) {
    return std::find(nodeMap.begin(),nodeMap.end(),node) != nodeMap.end();
}

class CalcPosB {
    State&   s;
    const LengthSet* lengthSet;
public:
    CalcPosB(State& ss, const LengthSet* lset)
        : s(ss), lengthSet(lset) {}
    Vector operator()(const Vector& pos) const
        { return lengthSet->calcPosB(s,pos); }
};

class CalcPosZ {
    const State&   s;
    const LengthSet* lengthSet;
public:
    CalcPosZ(const State& ss, const LengthSet* constraint) 
        : s(ss), lengthSet(constraint) {}
    Vector operator()(const Vector& b) const 
        { return lengthSet->calcPosZ(s, b); }
};

class CalcVelB {
    State&           s;
    const LengthSet* lengthSet;
public:
    CalcVelB(State& ss, const LengthSet* constraint)
        : s(ss), lengthSet(constraint) {}
    Vector operator()(const Vector& vel) 
        { return lengthSet->calcVelB(s,vel); }
};

//
// Calculate the position constraint violation (zero when constraint met).
//
Vector
LengthSet::calcPosB(State& s, const Vector& pos) const
{
    setPos(s, pos);

    // Although we're not changing the qErr cache entry here, we access
    // it with "upd" because setPos() will have modified the q's and
    // thus invalidated stage Position, so a "get" would fail.
    const Vector& qErr = getRBTree().updQErr(s);
    Vector b((int)loops.size());
    for (int i=0; i<(int)loops.size(); ++i)
        b[i] = loops[i].rbDistCons->getPosErr(qErr);
    return b;
}

//
// Calculate the velocity constraint violation (zero when constraint met).
//
Vector
LengthSet::calcVelB(State& s, const Vector& vel) const 
{
    setVel(s, vel);

    // Although we're not changing the uErr cache entry here, we access
    // it with "upd" because setVel() will have modified the u's and
    // thus invalidated stage Velocity, so a "get" would fail.
    const Vector& uErr = getRBTree().updUErr(s);

    Vector b((int)loops.size());
    for (int i=0; i<(int)loops.size(); ++i)
        b[i] = loops[i].rbDistCons->getVelErr(uErr);
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
// Let P = d perr/dq, A = d verr/du, Q = d qdot/du.
//
// We want to find a change deltaq that will elimate the current error b:
// P deltaq = b. [WRONG:]Instead we solve A * x = b, where x = inv(Q) * deltaq,
// and then solve deltaq = Q * x. Conveniently Q is our friendly invertible
// relation between qdot's and u's: qdot = Q*u.
//
// TODO: Sherm 080101: the above is incorrect. We have to factor (PQ^-1) and
// solve (PQ^-1)*deltaq = b for LS deltaq, which is not the same as 
// solving for LS x in P*x=b and then deltaq=Q*x.
//
// TODO: I have oversimplified the above since we are really looking for
// a least squares solution to an underdetermined system. With the 
// pseudoinverse A+ available we can write x = A+ * b.
//
Vector
LengthSet::calcPosZ(const State& s, const Vector& b) const
{
	const Matrix Gt = calcGrad(s);
    const Vector x = calcPseudoInverseA(Gt) * b;

    const SBModelVars&     mv = getRBTree().getModelVars(s);
    const Vector&          q  = getRBTree().getQ(s);
    const SBPositionCache& cc = getRBTree().updPositionCache(s);

    Vector       zu(getRBTree().getTotalDOF(),0.);
    Vector       zq(getRBTree().getTotalQAlloc(),0.);

    // map the vector dir back to the appropriate elements of z
    int indx=0; // sherm 060222: I added this
    for (int i=0 ; i<(int)nodeMap.size() ; i++) {
        const int d    = nodeMap[i]->getDOF();
        const int offs = nodeMap[i]->getUIndex();
        zu(offs,d) = x(indx,d);
        indx += d;

        // Make qdot = Q*u.
        nodeMap[i]->calcQDot(mv,q,cc,zu,zq);  // change u's to qdot's
    }
    assert(indx == ndofThisSet);

    return zq;
}

//
// Given a vector containing violations of velocity constraints, calculate
// a state update which would drive those violations to zero (if they
// are linear and well conditioned).
//
// This is simpler than CalcPosZ because we have the right gradient here (it's
// the same one in both places). TODO: and shouldn't be recalculated!
//
class CalcVelZ {
    const State&   s;
    const LengthSet* lengthSet;
	const Matrix Gt;
    const Matrix GInverse;
public:
    CalcVelZ(const State& ss, const LengthSet* lset)
      : s(ss), lengthSet(lset), Gt(lset->calcGrad(s)),
        GInverse(LengthSet::calcPseudoInverseA(Gt)) 
    {
    }

    Vector operator()(const Vector& b) {
        const Vector dir = GInverse * b;
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
LengthConstraints::enforcePositionConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
{
    assert(rbTree.getStage(s) >= Stage::Position-1);
    Vector& pos = rbTree.updQ(s);

    bool anyChanges = false;

    try { 
        for (int i=0 ; i<(int)pvConstraints.size() ; i++) {
            anyChanges = true; // TODO: assuming for now
            posMin.calc(requiredTol, desiredTol, pos,
                        CalcPosB(s, &pvConstraints[i]),
                        CalcPosZ(s, &pvConstraints[i]));
        }
    }
    catch ( SimTK::Exception::NewtonRaphsonFailure cptn ) {
        cout << "LengthConstraints::enforcePositionConstraints: exception: "
             << cptn.getMessage() << '\n';
    } 

    return anyChanges;
}

// Project out the velocity constraint errors from the given state. 
bool
LengthConstraints::enforceVelocityConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
{
    assert(rbTree.getStage(s) >= Stage(Stage::Velocity).prev());
    Vector& vel = rbTree.updU(s);

    bool anyChanges = false;

    try { 
        for (int i=0 ; i<(int)pvConstraints.size() ; i++) {
            anyChanges = true; // TODO: assuming for now
            velMin.calc(requiredTol, desiredTol, vel,
                        CalcVelB(s, &pvConstraints[i]),
                        CalcVelZ(s, &pvConstraints[i]));
        }
    }
    catch ( SimTK::Exception::NewtonRaphsonFailure cptn ) {
        cout << "LengthConstraints::enforceVelocityConstraints: exception: "
             << cptn.getMessage() << '\n';
    } 

    return anyChanges;
}

//
//// on each iteration
//// for each loop
//// -first calc all usual properties
//// -recursively compute phi_ni for each length constraint
////   from tip to outmost common body of 
//// -compute gradient
//// -update theta using quasi-Newton-Raphson
//// -compute Cartesian coords
//// -check for convergence
//// -repeat
//


void LengthSet::setPos(State& s, const Vector& pos) const
{
    const SBModelVars& mv = getRBTree().getModelVars(s);
    Vector&            q  = getRBTree().updQ(s);
    SBPositionCache&   pc = getRBTree().updPositionCache(s);
    Vector&            qErr = getRBTree().updQErr(s);

    for (int i=0 ; i<(int)nodeMap.size() ; i++)
        nodeMap[i]->copyQ(mv, pos, q);

    // TODO: sherm this is the wrong place for the stage update!
    s.invalidateAll(Stage::Position);

    // sherm TODO: this now computes kinematics for the whole system,
    // but it should only have to update the loop we are interested in.
    // Schwieters had this right before because his equivalent of 'setQ'
    // above also performed the kinematics, while ours just saves the
    // new state variable values and calculates here:
    getRBTree().realizeSubsystemPosition(s);

    // TODO: This is redundant after realizePosition(), but I'm leaving
    // it here because this is actually all that need be recalculated for
    // the loop closure iterations.
    for (int i=0; i<(int)loops.size(); ++i)
        loops[i].calcPosInfo(qErr,pc);
}

// Must have called LengthSet::setPos() already.
void LengthSet::setVel(State& s, const Vector& vel) const
{
    const SBModelVars&     mv = getRBTree().getModelVars(s);
    const SBPositionCache& pc = getRBTree().getPositionCache(s);
    Vector&                u  = getRBTree().updU(s);
    SBVelocityCache&       vc = getRBTree().updVelocityCache(s);
    Vector&                uErr = getRBTree().updUErr(s);

    for (int i=0 ; i<(int)nodeMap.size() ; i++)
        nodeMap[i]->copyU(mv, vel, u);

    // TODO: sherm this is the wrong place for the stage update!
    s.invalidateAll(Stage::Velocity);

    getRBTree().realizeSubsystemVelocity(s);

    // TODO: see comment above in setPos
    for (int i=0; i<(int)loops.size(); ++i)
        loops[i].calcVelInfo(pc,uErr,vc);
}

//
////A = df / dtheta_i
////  = [ -(q_ni-qn0)x , 1 ] [ phi_na^T Ha , phi_nb^T Hb , ... ]^T
////
//

// Calculate gradient by central difference for testing the analytic
// version. Presumes that calcEnergy has been called previously with current 
// value of ipos.
void 
LengthSet::fdgradf(State& s,
                   const Vector&  pos,
                   Matrix&        grad) const 
{
    const SBModelVars& mv = getRBTree().getModelVars(s);

    // Gradf gradf(tree);
    // gradf(x,grad); return;
    const double eps = 1e-6;
    const CalcPosB calcB(s, this);
    const Vector b = calcB(pos);
    int grad_indx=0;
    for (int i=0 ; i<(int)nodeMap.size() ; i++) {
        int pos_indx=nodeMap[i]->getQIndex();
        for (int j=0 ; j<nodeMap[i]->getNQInUse(mv) ; j++,pos_indx++,grad_indx++) {
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
LengthSet::testGrad(State& s, const Vector& pos, const Matrix& grad) const
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
// sherm 060222: OK, this routine doesn't really do what it says
// when the constraint includes a ball or free joint. It
// does not use the actual q's, which are quaternions. Instead
// it is something like d(v+ - v-)/du. If the u's are the 
// angular velocity across the joint then this is
// d(p+ - p-)/d qbar where qbar is a 1-2-3 Euler sequence
// which is 0,0,0 at the current orientation. ("instant
// coordinates").
//
// sherm: 060303 This routine uses the joint transition matrices ~H,
// which can be thought of as Jacobians ~H = d V_PB_G / d uB, that is
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
LengthSet::calcGrad(const State& s) const
{
    // We're not updating, but need to use upd here because Position stage
    // was invalidated by change to state.
    const SBPositionCache& pc = getRBTree().updPositionCache(s);

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
                    phiT[b][j] = phiT[b][j+1] * ~n->getPhi(pc);
                }
            }
        }

        // compute gradient
        //Vec3 uBond = unitVec(l.tipPos(pc,2) - l.tipPos(pc,1));
        const Vec3 uBond = l.tipPos(pc,2) - l.tipPos(pc,1); // p
        Row<2,Mat33> J[2];
        for (int b=1 ; b<=2 ; b++)
            // TODO: get rid of this b-1; make tips 0-based
            J[b-1] = Row<2,Mat33>(-crossMat(l.tipPos(pc,b) -
                                            l.tips(b).getNode().getX_GB(pc).T()),   one);
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
                const SpatialVec& HtCol = ~nodeMap[j]->getHRow(pc, k);
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
/*static*/ Matrix
LengthSet::calcPseudoInverseA(const Matrix& transposeOfA)
{
    const Matrix A = ~transposeOfA; // now A is dg/dtheta
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

   
//
// Given a vector in mobility space, project it along the motion constraints
// of this LengthSet, by removing its component normal to the constraint
// manifold. Here we want to obtain a least squares solution x to
//   A x = A v, A=[ncXnu] is the constraint Jacobian.
// We solve with pseudo inverse (see calcPseudoInverseA):
//   Least squares x = A+ * Av.
// That is the component of v which is normal to the constraint manifold.
// So we then set v = v - x and return.
//
// TODO: projections for integration errors should be done in the
// error norm. I'm not sure whether that has to be done in this routine
// or whether caller can scale on the way in and out.
//
void
LengthSet::projectUVecOntoMotionConstraints(const State& s, Vector& v)
{
    const Matrix transposeOfA = calcGrad(s);
    const Matrix pinvA = calcPseudoInverseA(transposeOfA); // TODO: this should already be available

    const Vector rhs  = packedMatTransposeTimesVec(transposeOfA,v); // i.e., rhs = Av

    //FIX: using inverse is inefficient
    const Vector x = pinvA * rhs;

    // subtract forces due to these constraints
    subtractPackedVecFromVec(v, x); 
}

Matrix
LengthSet::calcPseudoInverseAFD(const State& s) const
{
    Matrix grad(ndofThisSet,loops.size()); // <-- transpose of the actual dg/dtheta
    State sTmp = s;
    Vector pos = getRBTree().getQ(s); // a copy
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
//  1) calculate Y (block diagonal matrix) for all nodes (common body to tip)
//  2) for each LengthSet, calculate Lagrange multiplier(s) and add in
//     resulting force- update accelerations
//     Do this step from tip to base.
//

bool LengthConstraints::calcConstraintForces(const State& s, const Vector& udotErr,
                                             Vector& multipliers, 
                                             SBAccelerationCache& ac) const 
{
    if ( accConstraints.size() == 0 )
        return false;

    rbTree.calcY(s); // TODO <-- this doesn't belong here!

    for (int i=accConstraints.size()-1 ; i>=0 ; i--)
        accConstraints[i].calcConstraintForces(s,udotErr,multipliers,ac);

    return true;
}

void LengthConstraints::addInCorrectionForces(const State& s, const SBAccelerationCache& ac,
                                              SpatialVecList& spatialForces) const 
{
    for (int i=accConstraints.size()-1 ; i>=0 ; i--)
        accConstraints[i].addInCorrectionForces(s, ac, spatialForces);
}

void 
LengthConstraints::projectUVecOntoMotionConstraints(const State& s, Vector& vec)
{
    if ( pvConstraints.size() == 0 )
        return;

    for (int i=pvConstraints.size()-1 ; i>=0 ; i--)
        pvConstraints[i].projectUVecOntoMotionConstraints(s, vec);

    for (int i=pvConstraints.size()-1 ; i>=0 ; i--) 
        pvConstraints[i].testProjectedVec(s, vec);
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
static Real
computeA(const SBPositionCache& cc, 
         const SBDynamicsCache&      dc,
         const Vec3&    v1,
         const LoopWNodes& loop1, int s1,
         const LoopWNodes& loop2, int s2,
         const Vec3&    v2)
{
    const RigidBodyNode* n1 = &loop1.tips(s1).getNode();
    const RigidBodyNode* n2 = &loop2.tips(s2).getNode();

    const Mat33 one(1);

    SpatialRow t1 = ~v1 * Row<2,Mat33>(crossMat(n1->getX_GB(cc).T() - loop1.tipPos(cc,s1)), one);
    SpatialVec t2 = Vec<2,Mat33>(crossMat(loop2.tipPos(cc,s2) - n2->getX_GB(cc).T()), one) * v2;

    while ( n1->getLevel() > n2->getLevel() ) {
        t1 = t1 * ~n1->getPsi(dc);
        n1 = n1->getParent();
    }

    while ( n2->getLevel() > n1->getLevel() ) {
        t2 = n2->getPsi(dc) * t2;
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
        t1 = t1 * ~n1->getPsi(dc);
        t2 = n2->getPsi(dc) * t2;
        n1 = n1->getParent();
        n2 = n2->getParent();
    }

    // here n1==n2

    Real ret = t1 * n1->getY(dc) * t2;
    return ret;
}

//
// To be called for LengthSets consecutively from tip to base.
// Given acceleration errors for each loop contained in this LengthSet,
// this will calculate a force for every station, and store that force
// in the runtime block associated with that loop in the AccelerationCache
// output argument. It is up to the caller to do something with these forces.
//
// See Section 2.6 on p. 294 of Schwieters & Clore, 
// J. Magnetic Resonance 152:288-302. Equation reference below
// are to that paper. (sherm)
//
void
LengthSet::calcConstraintForces(const State& s, const Vector& udotErr,
                                Vector& multipliers, SBAccelerationCache& ac) const
{ 
    const SBPositionCache& pc      = getRBTree().getPositionCache(s);
    const SBVelocityCache& vc      = getRBTree().getVelocityCache(s);
    const SBDynamicsCache& dc      = getRBTree().getDynamicsCache(s);

    // This is the acceleration error for each loop constraint in this
    // LengthSet. We get a single scalar error per loop, since each
    // contains one distance constraint.
    // See Eq. [53] and the last term of Eq. [66].
    Vector rhs(loops.size());
    for (int i=0; i<(int)loops.size(); ++i)
        rhs[i] = loops[i].rbDistCons->getAccErr(udotErr);

    // Here A = Q*(J inv(M) J')*Q' where J is the kinematic Jacobian for
    // the constrained points and Q is the constraint Jacobian. See first 
    // term of Eq. [66].
    Matrix A(loops.size(),loops.size(),0.);
    for (int i=0 ; i<(int)loops.size() ; i++) {
        const Vec3 v1 = loops[i].tipPos(pc,2) - loops[i].tipPos(pc,1);
        for (int bi=1 ; bi<=2 ; bi++)
            for (int bj=1 ; bj<=2 ; bj++) {
                for (int j=i ; j<(int)loops.size() ; j++) {
                    const Vec3 v2 = loops[j].tipPos(pc,2) - loops[j].tipPos(pc,1);
                    Real  contrib = computeA(pc, dc, v1, loops[i], bi,
                                                     loops[j], bj, v2);
                    A(i,j) += contrib * (bi==bj ? 1 : -1);
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
    const Vector multipliersForThisSet = A.invert() * rhs;

    //TODO: need to copy out the right multipliers to the full multipliers array

    cout << "  OLD lambda = " << multipliersForThisSet << endl;

    // add forces due to these constraints
    for (int i=0 ; i<(int)loops.size() ; i++) {
        const Vec3 frc = multipliersForThisSet[i] * (loops[i].tipPos(pc,2) - loops[i].tipPos(pc,1));
        loops[i].setTipForce(ac, 2, -frc);
        loops[i].setTipForce(ac, 1,  frc);
    }
}

void LengthSet::addInCorrectionForces(const State& s, const SBAccelerationCache& ac,
                                      SpatialVecList& spatialForces) const 
{
    const SBPositionCache& pc = getRBTree().getPositionCache(s);

    for (int i=0; i<(int)loops.size(); ++i) {
        for (int t=1; t<=2; ++t) {
            const RigidBodyNode& node = loops[i].tipNode(t);
            const Vec3 force = loops[i].tipForce(ac,t);
            const Vec3 moment = cross(loops[i].tipPos(pc,t) - node.getX_GB(pc).T(), force);
            spatialForces[node.getNodeNum()] += SpatialVec(moment, force);
        }
    }
}

void LengthSet::testAccel(const State& s) const
{
    const Vector& udotErr = getRBTree().getUDotErr(s);

    double testTol=1e-8;
    for (int i=0 ; i<(int)loops.size() ; i++) {
        const Real aerr = fabs(loops[i].rbDistCons->getAccErr(udotErr));
        if (aerr > testTol)
            cout << "LengthSet::testAccel: constraint condition between atoms "
                 << loops[i].tips(1) << " and " << loops[i].tips(2) << " violated.\n"
                 << "\tnorm of violation: " << aerr << '\n';
    }
    cout.flush();
}


// This just computes ~mat*vec but first plucks out the relevant entries
// in vec to squash it down to the same size as mat.
Vector 
LengthSet::packedMatTransposeTimesVec(const Matrix& packedMat, const Vector& vec)
{
    // sherm 060222: 
    assert(vec.size() == getRBTree().getTotalDOF());
    assert(packedMat.nrow() == ndofThisSet && packedMat.ncol() == ndofThisSet);

    Vector packedVec(ndofThisSet);    //build vector same size as mat
    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        packedVec(indx,d) = vec(offs,d);
        indx += d;
    }
    // sherm 060222: 
    assert(indx == ndofThisSet);

    return ~packedMat * packedVec; 
}

// vec is a full vector in mobility space; packedVec consists of just
// the mobilities used in this LengthSet.
void
LengthSet::subtractPackedVecFromVec(Vector& vec,
                                    const Vector& packedVec)
{
    // sherm 060222:
    assert(vec.size() == getRBTree().getTotalDOF());
    assert(packedVec.size() == ndofThisSet);

    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        vec(offs,d) -= packedVec(indx,d);
        indx += d;
    }

    assert(indx == ndofThisSet);
}

void
LengthSet::testProjectedVec(const State& s, 
                            const Vector& vec) const
{
    // sherm 060222:
    assert(vec.size() == getRBTree().getTotalDOF());

    Vector packedVec(ndofThisSet);
    int indx=0;
    for (int j=0 ; j<(int)nodeMap.size() ; j++) {
        const int d    = nodeMap[j]->getDOF();
        const int offs = nodeMap[j]->getUIndex();
        packedVec(indx,d) = vec(offs,d);
        indx += d;
    }
    assert(indx == ndofThisSet);

    const Matrix grad = calcGrad(s);
    const Vector test = ~grad * packedVec;

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
LengthConstraints::fixVel0(State& s, Vector& iVel)
{
    // use molecule grouping
    for (int i=0 ; i<(int)accConstraints.size() ; i++)
        accConstraints[i].fixVel0(s, iVel);
}

//
// Correct internal velocities such that the constraint conditions are
// met under the condition that the disturbance to station velocites is
// small as possible.
//
void
LengthSet::fixVel0(State& s, Vector& iVel)
{
    assert(iVel.size() == getRBTree().getTotalDOF());

    const SBPositionCache& pc = getRBTree().getPositionCache(s);
    const Vector&          uErr = getRBTree().getUErr(s);

    // store internal velocities
    Vector iVel0 = iVel;

    // Allocate the vector of body impulses so we don't have to do it
    // inside the loop below. No need to initialize them here, though.
    Vector_<SpatialVec> bodyImpulses(getRBTree().getNBodies());

    //TODO
    // Currently we do not have any constraints involving mobilities (u's)
    // directly, so we don't generate mobility impulses here. But we'll need
    // a set of zero mobility forces to apply below.
    Vector zeroMobilityImpulses(getRBTree().getTotalDOF());
    zeroMobilityImpulses.setToZero();

    // verr stores the current velocity errors, which we're assuming are valid.
    Vector verr((int)loops.size());
    for (int i=0; i<(int)loops.size(); ++i)
        verr[i] = loops[i].rbDistCons->getVelErr(uErr);

    Matrix mat(loops.size(),loops.size());
    std::vector<Vector> deltaIVel(loops.size());
    for (int m=0 ; m<(int)loops.size() ; m++) {
        deltaIVel[m].resize(iVel.size());

        // Set all velocities to zero. TODO: this should just be an "ignore velocity"
        // option to realizeVelocity(); it shouldn't actually require putting zeroes everywhere.
        iVel = 0.;
        getRBTree().setU(s, iVel );
        getRBTree().realizeSubsystemVelocity(s);

        // sherm: I think the following is a unit "probe" velocity, projected
        // along the separation vector. 
        // That would explain the fact that there are no velocities here!
        const Vec3 probeImpulse = loops[m].tipPos(pc,2)-loops[m].tipPos(pc,1);
        const Vec3 force1 = -probeImpulse;
        const Vec3 force2 =  probeImpulse;

        const RigidBodyNode& node1 = loops[m].tipNode(1);
        const RigidBodyNode& node2 = loops[m].tipNode(2);
        const Vec3 moment1 = cross(loops[m].tipPos(pc,1)-node1.getX_GB(pc).T(), force1);
        const Vec3 moment2 = cross(loops[m].tipPos(pc,2)-node2.getX_GB(pc).T(), force2);

        // Convert the probe impulses at the stations to spatial impulses at
        // the body origin.
        bodyImpulses.setToZero();
        bodyImpulses[node1.getNodeNum()] += SpatialVec(moment1, force1);
        bodyImpulses[node2.getNodeNum()] += SpatialVec(moment2, force2);

        // Apply probe impulse as though it were a force; the resulting "acceleration"
        // is actually the deltaV produced by this impulse, that is, a deltaV which
        // results in a change in velocity along the line between the stations.

        // calc deltaVa from k-th constraint condition
        //FIX: make sure that nodeMap is ordered by level
        // get all nodes in given molecule, ordered by level
        getRBTree().calcZ(s, zeroMobilityImpulses, bodyImpulses);
        getRBTree().calcTreeAccel(s);
        deltaIVel[m] = getRBTree().getUDot(s);

        // Set the new velocities.
        getRBTree().setU(s, deltaIVel[m] );
        getRBTree().realizeSubsystemVelocity(s);

        // Calculating partial(velocityError[n])/partial(deltav[m]). Any velocity
        // we see here is due to the deltav, since we started out at zero (uErr
        // was modified by realizeVelocity(), but our reference is still valid).
        for (int n=0; n<(int)loops.size(); ++n)
            mat(n,m) = loops[n].rbDistCons->getVelErr(uErr);

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
    getRBTree().realizeSubsystemVelocity(s);
}


std::ostream& operator<<(std::ostream& o, const RBStation& s) {
    o << "station " << s.getPoint() << " on node " << s.getNode().getNodeNum();
    return o;
}


std::ostream& operator<<(std::ostream& o, const RBDirection& d) {
    o << "normal " << d.getUnitVec() << " on node " << d.getNode().getNodeNum();
    return o;
}

    ////////////////////////////
    // RB DISTANCE CONSTRAINT //
    ////////////////////////////

void RBDistanceConstraint::calcStationPosInfo(int i, 
        SBPositionCache& pc) const
{
    updStation_G(pc,i) = getNode(i).getX_GB(pc).R() * getPoint(i);
    updPos_G(pc,i)     = getNode(i).getX_GB(pc).T() + getStation_G(pc,i);
}

void RBDistanceConstraint::calcStationVelInfo(int i, 
        const SBPositionCache& pc, 
        SBVelocityCache&       vc) const
{
    const Vec3& w_G = getNode(i).getSpatialAngVel(vc);
    const Vec3& v_G = getNode(i).getSpatialLinVel(vc);
    updStationVel_G(vc,i) = cross(w_G, getStation_G(pc,i));
    updVel_G(vc,i)        = v_G + getStationVel_G(vc,i);
}

void RBDistanceConstraint::calcStationAccInfo(int i, 
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        SBAccelerationCache&   ac) const
{
    const Vec3& w_G  = getNode(i).getSpatialAngVel(vc);
    const Vec3& v_G  = getNode(i).getSpatialLinVel(vc);
    const Vec3& aa_G = getNode(i).getSpatialAngAcc(ac);
    const Vec3& a_G  = getNode(i).getSpatialLinAcc(ac);
    updAcc_G(ac,i) = a_G + cross(aa_G, getStation_G(pc,i))
                         + cross(w_G,  getStationVel_G(vc,i)); // i.e., w X (wXr)
}


void RBDistanceConstraint::calcPosInfo(Vector& qErr, SBPositionCache& pc) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationPosInfo(i,pc);

    const Vec3 p = getPos_G(pc,2) - getPos_G(pc,1);
    updFromTip1ToTip2_G(pc)  = p;
    const Real separation  = getFromTip1ToTip2_G(pc).norm();
    updUnitDirection_G(pc)   = getFromTip1ToTip2_G(pc) / separation;
    //TODO:  |p|-d (should be 0.5(p^2-d^2)
    //updPosErr(cc) = separation - distance; 
    updPosErr(qErr) = 0.5*(p.normSqr() - distance*distance);

}

void RBDistanceConstraint::calcVelInfo(
        const SBPositionCache& pc, 
        Vector&                uErr,
        SBVelocityCache&       vc) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationVelInfo(i,pc,vc);

    updRelVel_G(vc) = getVel_G(vc,2) - getVel_G(vc,1);
    //TODO: u.v, u=p/|p| (should be p.v)
    //updVelErr(mc)   = ~getUnitDirection_G(pc) * getRelVel_G(vc);
    updVelErr(uErr) = ~getFromTip1ToTip2_G(pc) * getRelVel_G(vc);
}

void RBDistanceConstraint::calcAccInfo(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        Vector&                udotErr,
        SBAccelerationCache&   ac) const
{
    assert(isValid() && distConstNum >= 0);
    for (int i=1; i<=2; ++i) calcStationAccInfo(i,pc,vc,ac);

    // TODO: v.v + a.p (good), but would have to be
    // u.a + (v-(v.u).u).v/|p| to be compatible with above
    const Vec3 relAcc_G = getAcc_G(ac,2) - getAcc_G(ac,1);
    updAccErr(udotErr) = getRelVel_G(vc).normSqr() + (~relAcc_G * getFromTip1ToTip2_G(pc));
}

