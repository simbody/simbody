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

#include "IVMLengthConstraints.h"
#include "IVMRigidBodyNode.h"
#include "IVMRigidBodyTree.h"
#include "cdsVec3.h"
#include "cdsListAutoPtr.h"
#include "RVec.h"
#include "cdsMath.h"
#include "cdsMatrix.h"
#include "fixedVector.h"
#include "cdsVector.h"
#include "subVector.h"
#include "subMatrix.h"
#include "matrixTools.h"
#include "cdsNewtonRaphson.h"
#include "cdsIomanip.h"

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using MatrixTools::inverse;
#endif /* USE_CDS_NAMESPACE */

typedef CDSMatrix<double>     RMat;
typedef CDSVector<double,0>   RVec0;
typedef SubVector<RVec>       SubVec;
typedef SubVector<const RVec> ConstSubVec;
typedef SubMatrix<RMat>       SubMat;

class IVMLoopWNodes;
static int compareLevel(const IVMLoopWNodes& l1,    //forward declarations
                        const IVMLoopWNodes& l2);

//static bool sameBranch(const IVMRigidBodyNode* tip,
//                       const IVMLoopWNodes& l );

class BadNodeDef {};  //exception

/**
 * Collect up useful information about a loop. 
 * This includes the two connected atoms, ordered by level, and the
 * paths from each of the associated nodes back to the common ancestor.
 * We also identify the molecule base node for the molecule which
 * contains both ends of the loop.
 * We will throw an exception if the loop ends are both on the same
 * node or if they are on different molecules.
 */
class IVMLoopWNodes {
public:
    IVMLoopWNodes() : rbDistCons(0), rt(0), flipStations(false), base(0), moleculeNode(0) {}
    IVMLoopWNodes(const IVMDistanceConstraint&, CDSList<IVMDistanceConstraintRuntime>&);

    void calcPosInfo() const { rbDistCons->calcPosInfo(*rt); }
    void calcVelInfo() const { rbDistCons->calcVelInfo(*rt); }
    void calcAccInfo() const { rbDistCons->calcAccInfo(*rt); }

    const double& getDistance() const { return rbDistCons->getDistance(); }

    // Return one of the stations, ordered such that tips(1).level <= tips(2).level.
    const IVMStation& tips(int i) const {return rbDistCons->getStation(ix(i));}
    const IVMRigidBodyNode& tipNode(int i) const {return tips(i).getNode();}

    const CDSVec3& tipPos(int i)   const {return getTipRuntime(i).pos_G;}
    const CDSVec3& tipVel(int i)   const {return getTipRuntime(i).vel_G;}
    const CDSVec3& tipAcc(int i)   const {return getTipRuntime(i).acc_G;}
    const CDSVec3& tipForce(int i) const {return getTipRuntime(i).force_G;}

    // Use this for both forces and impulses.
    void setTipForce(int i, const CDSVec3& f) const { updTipRuntime(i).force_G = f; }

private:
    int ix(int i) const { assert(i==1||i==2); return flipStations ? 3-i : i; }
    const IVMStationRuntime& getTipRuntime(int i) const
      { return rt->stationRuntimes[ix(i)-1]; }
    IVMStationRuntime& updTipRuntime(int i) const
      { return rt->stationRuntimes[ix(i)-1]; }

    const IVMDistanceConstraint*  rbDistCons;  // reference to the constraint
    IVMDistanceConstraintRuntime* rt;          // ... and its runtime

    // calculated info about the constraint
    bool                           flipStations; // make sure station(1).level
                                                 //   <= station(2).level
    FixedVector<IVMNodePtrList,2,1> nodes;   // the two paths: base..tip1, base..tip2,
                                            //   incl. tip nodes but not base
    IVMRigidBodyNode*                 base;    // highest-level common ancestor of tips
    const IVMRigidBodyNode*           moleculeNode;

    friend class IVMLengthSet;
    friend ostream& operator<<(ostream& os, const IVMLengthSet& s);
    friend void IVMLengthConstraints::construct(CDSList<IVMDistanceConstraint>&,
                                                CDSList<IVMDistanceConstraintRuntime>&);
    friend int compareLevel(const IVMLoopWNodes& l1,
                            const IVMLoopWNodes& l2);
    friend bool sameBranch(const IVMRigidBodyNode* tip,
                           const IVMLoopWNodes& l );
};

IVMLoopWNodes::IVMLoopWNodes(const IVMDistanceConstraint& dc,
                       CDSList<IVMDistanceConstraintRuntime>& rts)
  : rbDistCons(&dc), rt(&rts[dc.getRuntimeIndex()]),flipStations(false), base(0), moleculeNode(0)
{
    using IVMInternalDynamics::IVMException;

    const IVMRigidBodyNode* dcNode1 = &dc.getStation(1).getNode();
    const IVMRigidBodyNode* dcNode2 = &dc.getStation(2).getNode();

    if (dcNode1==dcNode2) {
        cout << "IVMLoopWNodes::IVMLoopWNodes: bad topology:\n\t"
             << "loop stations " << dc.getStation(1)
             << " and  "         << dc.getStation(2)
             << " are now in the same node. Deleting loop.\n";
        throw BadNodeDef();
    }

    // Ensure that tips(2) is the atom which is farther from the base.
    flipStations = (dcNode1->getLevel() > dcNode2->getLevel());

    // OK to use tips() at this point.

    // Collect up the node path from tips(2) down to the last node on its
    // side of the loop which is at a higher level than tips(1) (may be none).
    IVMRigidBodyNode* node1 = &tips(1).getNode();
    IVMRigidBodyNode* node2 = &tips(2).getNode();
    while ( node2->getLevel() > node1->getLevel() ) {
        nodes(2).append(node2);
        node2 = node2->getParent();
    }

    // We're at the same level on both sides of the loop. Run down both
    // simultaneously until we hit the first common ancestor, collecting
    // up the nodes along the two paths (but not the common ancestor).
    while ( node1 != node2 ) {
        if ( node1->isGroundNode() ) {
            cerr << "IVMLoopWNodes::IVMLoopWNodes: could not find base node.\n\t"
                 << "loop between stations " << tips(1) << " and " 
                 << tips(2) << "\n";
            throw IVMException("IVMLoopWNodes::IVMLoopWNodes: could not find base node");
        }
        nodes(1).append(node1);
        nodes(2).append(node2);
        node1 = node1->getParent();
        node2 = node2->getParent();
    }

    base = node1;   // that's the common ancestor

    // We want these in base-to-tip order.
    nodes(1).reverse();
    nodes(2).reverse();

    // find molecule node (level==1 node)
    for (moleculeNode=base; 
         moleculeNode->getLevel()>1 ; 
         moleculeNode=moleculeNode->getParent()) ;

    if ( moleculeNode->getLevel()<1 ) {
        cerr << "IVMLoopWNodes::IVMLoopWNodes: could not find molecule node.\n\t"
             << "loop between atoms " << tips(1) << " and " 
             << tips(2) << "\n";
        throw IVMException("IVMLoopWNodes::IVMLoopWNodes: could not find molecule node");
    }
}

typedef CDSList<IVMLoopWNodes> IVMLoopList;

class IVMLengthSet {
    static void construct(const IVMLoopList& loops);
    const IVMLengthConstraints*  lConstraints;
    IVMRigidBodyTree&            rbTree;
    const int                    verbose;
    IVMLoopList                  loops;    
    int                          dim;
    CDSList<IVMRigidBodyNode*>   nodeMap; //unique nodes (union of loops->nodes)
public:
    IVMLengthSet(const IVMLengthConstraints* lConstraints, 
                 const IVMLoopWNodes& loop)
        : lConstraints(lConstraints), rbTree(lConstraints->rbTree), 
          verbose(lConstraints->verbose), dim(0)
    {
        addConstraint(loop);
    }

    void addConstraint(const IVMLoopWNodes& loop) {
        loops.append( IVMLoopWNodes(loop) );
        for (int b=1 ; b<=2 ; b++)
            for (int i=0; i<loop.nodes(b).size(); i++)
                if (!nodeMap.contains(loop.nodes(b)[i])) {
                    dim += loop.nodes(b)[i]->getDim();
                    nodeMap.append( loop.nodes(b)[i] );
                }
    }

    void determineCouplings();
    bool contains(IVMRigidBodyNode* node) {
        bool found=false;
        for (int i=0 ; i<loops.size() ; i++)
            if (    loops[i].nodes(1).contains(node)
                 || loops[i].nodes(2).contains(node) ) 
            {
                found=true;
                break;
            }
        return found;
    }

    void  setPos(const RVec& pos) const;
    void  setVel(const RVec& vel) const;
    RVec  getPos();
    RVec0 calcPosB(const RVec& pos) const;
    RVec0 calcVelB(const RVec& pos, const RVec& vel) const;
    RVec  calcPosZ(const RVec& b) const;
    RMat  calcGrad() const;
    RMat  calcGInverse() const;

    void  calcConstraintForces() const; // updates runtime only
    void  addInCorrectionForces(CDSVecVec6& spatialForces) const; // spatialForces+=correction

    void  fixVel0(RVec&);
    void  fixInternalForce(RVec&);
    RVec0 multiForce(const RVec&, const RMat& mat);
    void  addForce(RVec&, const RVec& ve);

    void testAccel();
    void testInternalForce(const RVec&);
    friend ostream& operator<<(ostream& os, const IVMLengthSet& s);
    friend class IVMCalcPosB;
    friend class IVMCalcVelB;
    friend class IVMCalcPosZ;
    friend class IVMCalcVelZ;

    void fdgradf(const RVec& pos, RMat& grad) const;
    void testGrad(const RVec& pos, const RMat& grad) const;
};

ostream& 
operator<<(ostream& os, const IVMLengthSet& s) 
{
    for (int i=0 ; i<s.loops.size() ; i++)
        os << setw(4) << s.loops[i].tips(1) << "<->" 
           << setw(4) << s.loops[i].tips(2)  << ": " 
           << s.loops[i].getDistance() << "  ";
    return os;
}

class IVMLengthConstraintsPrivates {
public:
    IVMLengthConstraintsPrivates() : posMin(cout), velMin(cout) {}

public:
    CDSListAutoPtr<IVMLengthSet> constraints;     // used for pos, vel
    CDSListAutoPtr<IVMLengthSet> accConstraints;  // used for acc
    CDSNewtonRaphson             posMin, velMin;
};

IVMLengthConstraints::IVMLengthConstraints(IVMRigidBodyTree& rbt, const double& ctol, int vbose)
    : bandCut(1e-7), maxIters( 20 ), maxMin( 20 ), rbTree(rbt), verbose(vbose)
{
    priv = new IVMLengthConstraintsPrivates;
    priv->posMin.maxIters = maxIters;
    priv->posMin.maxMin   = maxMin;
    priv->posMin.tol      = ctol;
    priv->velMin.maxIters = maxIters;
    priv->velMin.maxMin   = maxMin;
    priv->velMin.tol      = ctol*10;
}

IVMLengthConstraints::~IVMLengthConstraints()
{
    delete priv;
}

//
// compare IVMLoopWNodes by base node level value
//
static int
compareLevel(const IVMLoopWNodes& l1,
             const IVMLoopWNodes& l2) 
{ 
    if ( l1.base->getLevel() > l2.base->getLevel() ) 
        return 1;
    else if ( l1.base->getLevel() < l2.base->getLevel() )
        return -1;
    else 
        return 0;
}

//1) construct: given list of loops 
//   a) if appropriate issue warning and exit.
//   b) sort by base node of each loop.
//   c) find loops which intersect: combine loops and increment
//    number of length constraints
void
IVMLengthConstraints::construct(CDSList<IVMDistanceConstraint>& iloops,
                                CDSList<IVMDistanceConstraintRuntime>& rts)
{
    //clean up
    priv->constraints.resize(0);
    priv->accConstraints.resize(0);

    // if ( !useLengthConstraint )  FIX:!!!
    //   //issue error message?
    //   return;

    IVMLoopList loops;
    for (int i=0 ; i<iloops.size() ; i++) {
        try {
            IVMLoopWNodes loop( iloops[i], rts );
            loops.append( loop );
        }
        catch ( BadNodeDef ) {}
    }

    // sort loops by base->level
    loops.sort(compareLevel);
    IVMLoopList accLoops = loops;  //version for acceleration

    // find intersections -- this version keeps hierarchical loops distinct
    // sherm: the loops are considered coupled if a lower one includes the
    // first body up from the base along either branch of the upper one (if
    // it doesn't include that it can't include any further up the branch either).
    // This makes good sense to me.
    for (int i=0 ; i<loops.size() ; i++) {
        priv->constraints.append( new IVMLengthSet(this, loops[i]) );
        for (int j=i+1 ; j<loops.size() ; j++)
            for (int k=0 ; k<priv->constraints.size() ; k++)
                if (   (loops[j].nodes(1).size()
                        && priv->constraints[k]->contains(loops[j].nodes(1)[0]))
                    || (loops[j].nodes(2).size()
                        && priv->constraints[k]->contains(loops[j].nodes(2)[0])))
                {
                    //add length constraint to loop i
                    priv->constraints[i]->addConstraint(loops[j]);
                    loops.remove(j);
                    j--;
                    break;
                }
    }

    // find intersections - group all loops with tip->trunk relationship
    // TODO sherm: I don't understand why these have to be more coupled than
    // the pos/vel constraints. (This code couples all the loops on the same
    // molecule).
    for (int i=0 ; i<accLoops.size() ; i++) {
        priv->accConstraints.append( new IVMLengthSet(this, accLoops[i]) );
        for (int j=i+1 ; j<accLoops.size() ; j++)
            if ( accLoops[i].moleculeNode == accLoops[j].moleculeNode ) {
                priv->accConstraints[i]->addConstraint(accLoops[j]);
                accLoops.remove(j);
                j--;
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

    for (l_int i=0 ; i<priv->accConstraints.size() ; i++)
        priv->accConstraints[i]->determineCouplings(); // XXX not implemented

    if (priv->constraints.size()>0 && verbose&IVMInternalDynamics::printLoopInfo) {
        cout << "IVMLengthConstraints::construct: pos/vel length constraints found:\n";
        for (l_int i=0 ; i<priv->constraints.size() ; i++)
            cout << "\t" << *priv->constraints[i] << "\n";
        cout << "IVMLengthConstraints::construct: accel length constraints found:\n";
        for (l_int i=0 ; i<priv->accConstraints.size() ; i++)
            cout << "\t" << *priv->accConstraints[i] << "\n";
    }
}

class IVMCalcPosB {
    const IVMLengthSet* constraint;
public:
    IVMCalcPosB(const IVMLengthSet* constraint)
        : constraint(constraint) {}
    RVec0 operator()(const RVec& pos) const
        { return constraint->calcPosB(pos); }
};

class IVMCalcPosZ {
    const IVMLengthSet* constraint;
public:
    IVMCalcPosZ(const IVMLengthSet* constraint) : constraint(constraint) {}
    RVec0 operator()(const RVec& b) const 
        { return constraint->calcPosZ(b); }
};

class IVMCalcVelB {
    const IVMLengthSet* constraint;
    const RVec& pos;
public:
    IVMCalcVelB(const IVMLengthSet* constraint, const RVec& pos)
        : constraint(constraint), pos(pos) {}
    RVec0 operator()(const RVec& vel) 
        { return constraint->calcVelB(pos,vel); }
};

//
// calculate the constraint violation (zero when constraint met)
//
RVec0
IVMLengthSet::calcPosB(const RVec& pos) const
{
    setPos(pos);

    RVec0 b( loops.size() );
    for (int i=0 ; i<loops.size() ; i++) 
        b(i) = loops[i].getDistance() - 
                sqrt(abs2(loops[i].tipPos(1) - loops[i].tipPos(2)));
    return b;
}

//
// calculate the constraint violation (zero when constraint met)
//
RVec0
IVMLengthSet::calcVelB(const RVec& pos, const RVec& vel) const 
{
    setPos(pos); // TODO (sherm: this is probably redundant)
    setVel(vel);

    RVec0 b( loops.size() );
    for (int i=0 ; i<loops.size() ; i++) {
        // TODO why the minus sign here? (doesn't work right without it) sherm
        b(i) = -dot( unitVec(loops[i].tipPos(2) - loops[i].tipPos(1)),
                    loops[i].tipVel(2) - loops[i].tipVel(1) );
    }

    return b;
}

//
// Given a vector containing violations of position constraints, calculate
// a state update which should drive those violations to zero (if they
// were linear).
//
RVec
IVMLengthSet::calcPosZ(const RVec& b) const
{
    RVec dir = calcGInverse() * b;
    RVec z(rbTree.getDim(),0.0);

    // map the vector dir back to the appropriate elements of z
    for (int i=0 ; i<nodeMap.size() ; i++) 
        SubVec(z,nodeMap[i]->getStateOffset(), nodeMap[i]->getDim())
            = SubVec(dir,i+1,nodeMap[i]->getDim());

    return z;
}

//
// Given a vector containing violations of velocity constraints, calculate
// a state update which would drive those violations to zero (if they
// are linear and well conditioned).
//
class IVMCalcVelZ {
    const IVMLengthSet* constraint;
    const RMat       gInverse;
public:
    IVMCalcVelZ(const IVMLengthSet* constraint)
      : constraint(constraint), gInverse(constraint->calcGInverse()) {}
    RVec operator()(const RVec& b) {
        RVec dir = gInverse * b;
        RVec z(constraint->rbTree.getDim(),0.0);

        // map the vector dir back to the appropriate elements of z
        for (int i=0 ; i<constraint->nodeMap.size() ; i++) {
            const IVMRigidBodyNode* n = constraint->nodeMap[i];
            SubVec(z,n->getStateOffset(),n->getDim()) =  SubVec(dir,i+1,n->getDim());
        }
        return z;
    }
};

//
// Project out the position and velocity constraint errors from the given
// state. 
// XXX (sherm): Velocity errors should be calculated using the updated positions
// (that is, after the position correction. I don't think that is happening
// below.
void
IVMLengthConstraints::enforce(RVec& pos, RVec& vel)
{
    priv->posMin.verbose  = ((verbose&IVMInternalDynamics::printLoopDebug) != 0);
    priv->velMin.verbose  = ((verbose&IVMInternalDynamics::printLoopDebug) != 0);
    try { 
        for (int i=0 ; i<priv->constraints.size() ; i++) {
            if ( verbose&IVMInternalDynamics::printLoopDebug )
                cout << "IVMLengthConstraints::enforce: position " 
                     << *(priv->constraints[i]) << '\n';
            priv->posMin.calc(pos,
                              IVMCalcPosB(priv->constraints[i].get()),
                              IVMCalcPosZ(priv->constraints[i].get()));
        }
        for (int i=0 ; i<priv->constraints.size() ; i++) {
            if ( verbose&IVMInternalDynamics::printLoopDebug )
                cout << "IVMLengthConstraints::enforce: velocity " 
                     << *priv->constraints[i] << '\n';
            priv->velMin.calc(vel,
                              IVMCalcVelB(priv->constraints[i].get(),pos),
                              IVMCalcVelZ(priv->constraints[i].get()));
        }
    }
    catch ( CDSNewtonRaphson::Fail cptn ) {
        cout << "IVMLengthConstraints::enforce: exception: "
             << cptn.mess << '\n';
    } 
// catch ( ... ) {
//   cout << "IVMLengthConstraints::enforce: exception: "
//      << "uncaught exception in NewtonRaphson.\n" << ends;
//   cout.flush();
//   throw;
// }
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


void IVMLengthSet::setPos(const RVec& pos) const
{
    for (int i=0 ; i<nodeMap.size() ; i++)
        nodeMap[i]->setPos(pos); //also calc necessary properties

    for (int i=0; i<loops.size(); ++i)
        loops[i].calcPosInfo();
}

// Must have called IVMLengthSet::setPos() already.
void IVMLengthSet::setVel(const RVec& vel) const
{
    for (int i=0 ; i<nodeMap.size() ; i++)
        nodeMap[i]->setVel(vel); //also calc necessary properties

    for (int i=0; i<loops.size(); ++i)
        loops[i].calcVelInfo();
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
IVMLengthSet::fdgradf( const RVec& pos,
                    RMat&       grad) const 
{
    // Gradf gradf(tree);
    // gradf(x,grad); return;
    const double eps = 1e-8;
    const IVMCalcPosB calcB(this);
    const RVec0 b = calcB(pos);
    int grad_indx=0;
    for (int i=0 ; i<nodeMap.size() ; i++) {
        int pos_indx=nodeMap[i]->getStateOffset();
        for (int j=0 ; j<nodeMap[i]->getDim() ; j++,pos_indx++,grad_indx++) {
            RVec posp = pos;
            posp(pos_indx) += eps;
            const RVec0 bp = calcB(posp);
            for (int k=0 ; k<b.size() ; k++)
                grad(grad_indx,k) = -(bp(k)-b(k)) / eps;
        }
    }
}

void 
IVMLengthSet::testGrad(const RVec& pos, const RMat& grad) const
{
    double tol = 1e-4;
    // Costf costf(tree);
    // double f = costf(x);
    RMat fdgrad(dim,loops.size());
    fdgradf(pos,fdgrad);

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
// unitVec(q+ - q-) * d (q+ - q-) / d (theta_i)
// d g / d theta for all hingenodes in nodemap
//
// sherm: this appears to calculate the transpose of G
//
RMat
IVMLengthSet::calcGrad() const
{
    RMat grad(dim,loops.size(),0.0);
    CDSMat33 one(0.0); one.setDiag(1.0);  //FIX: should be done once

    for (int i=0 ; i<loops.size() ; i++) {
        const IVMLoopWNodes& l = loops[i];
        FixedVector<CDSList<CDSMat66>,2,1> phiT;
        for (int b=1 ; b<=2 ; b++) {
            phiT(b).resize( l.nodes(b).size() );
            if ( l.nodes(b).size() ) {
                phiT(b)[phiT(b).size()-1].set(0.0);
                phiT(b)[phiT(b).size()-1].setDiag(1.0);
                for (int j=l.nodes(b).size()-2 ; j>=0 ; j-- ) {
                    IVMRigidBodyNode* n = l.nodes(b)[j+1];
                    phiT(b)[j] = phiT(b)[j+1] * transpose( n->getPhi() );
                }
            }
        }

        // compute gradient
        CDSVec3 uBond = unitVec(l.tipPos(2) - l.tipPos(1));
        FixedVector<FixedMatrix<double,3,6>,2,1> J;
        for (int b=1 ; b<=2 ; b++)
            J(b) = blockMat12(-crossMat(l.tipPos(b) -
                                        l.tips(b).getNode().getOB_G()), one);
        int g_indx=0;
        for (int j=0 ; j<nodeMap.size() ; j++) {
            RMat H = MatrixTools::transpose(nodeMap[j]->getH());
            double elem=0.0;
            int l1_indx = l.nodes(1).getIndex(nodeMap[j]);
            int l2_indx = l.nodes(2).getIndex(nodeMap[j]);
            for (int k=0 ; k<H.ncol() ; k++) {
                CDSVec6 Hcol = subCol(H,k,0,5).vector();
                if ( l1_indx >= 0 ) { 
                    elem = -dot(uBond , CDSVec3(J(1) * phiT(1)[l1_indx]*Hcol));
                } else if ( l2_indx >= 0 ) { 
                    elem = dot(uBond , CDSVec3(J(2) * phiT(2)[l2_indx]*Hcol));
                }
                grad(g_indx++,i) = elem;
            }
        }
    }
    return grad;
} 

static double 
abs2(const RMat& m)
{
    double ret=0.;
    for (int i=0 ; i<m.nrow() ; i++)
        for (int j=0 ; j<m.ncol() ; j++)
            ret += CDSMath::sq( m(i,j) );
    return ret;
}

//
// Calculate generalized inverse which minimizes changes in soln vector.
// TODO (sherm) I THINK this is trying to create a pseudoinverse
// using normal equations which is numerically bad. Should use an SVD
// or QTZ factorization instead.
//     Want least squares solution x to G'x=b. Normal equations are
//     GG'x=Gb, x=inv(GG')Gb ??? [returning G inv(G'G) below] ???
//     ??? either this is wrong or I don't understand
RMat
IVMLengthSet::calcGInverse() const
{
    RMat grad = calcGrad(); // <-- appears to be transpose of the actual dg/dtheta
    if ( verbose & IVMInternalDynamics::printLoopDebug ) {
        RVec pos(rbTree.getDim()); rbTree.getPos(pos);
        testGrad(pos,grad);
    }

    RMat ret(grad.nrow(),grad.nrow(),0.0); // <-- wrong dimension ??? sherm TODO
    if ( abs2(grad) > 1e-10 ) 
        ret = grad * inverse( MatrixTools::transpose(grad)*grad );
    return ret;
}

//acceleration:
//  0) after initial acceleration calculation:
//  1) calculate Y (block diagonal matrix) for all nodes (base to tip)
//  2) for each IVMLengthSet, calculate Lagrange multiplier(s) and add in
//     resulting force- update accelerations
//     Do this step from tip to base.
//

bool IVMLengthConstraints::calcConstraintForces() const {
    if ( priv->accConstraints.size() == 0 )
        return false;

    rbTree.calcY();

    for (int i=priv->accConstraints.size()-1 ; i>=0 ; i--)
        priv->accConstraints[i]->calcConstraintForces();

    return true;
}

void IVMLengthConstraints::addInCorrectionForces(CDSVecVec6& spatialForces) const {
    for (int i=priv->accConstraints.size()-1 ; i>=0 ; i--)
        priv->accConstraints[i]->addInCorrectionForces(spatialForces);
}

void 
IVMLengthConstraints::fixGradient(RVec& forceInternal)
{
    if ( priv->constraints.size() == 0 )
        return;

    for (int i=priv->constraints.size()-1 ; i>=0 ; i--)
        priv->constraints[i]->fixInternalForce(forceInternal);

    for (l_int i=priv->constraints.size()-1 ; i>=0 ; i--) 
        priv->constraints[i]->testInternalForce(forceInternal);
}

typedef SubVector<const CDSVec6>       SubVec6;
typedef SubVector<const CDSVec6>       SubVec6;

//
// Calculate the acceleration of atom assuming that the spatial acceleration
// of the body (node) it's on is available.
//
//static CDSVec3
//getAccel(const IVMAtom* a)
//{
//    const IVMRigidBodyNode* n = a->node;
//    CDSVec3 ret( SubVec6(n->getSpatialAcc(),3,3).vector() );
//    ret += cross( SubVec6(n->getSpatialAcc(),0,3).vector() , a->pos - n->getAtom(0)->pos );
//    ret += cross( SubVec6(n->getSpatialVel(),0,3).vector() , a->vel - n->getAtom(0)->vel );
//    return ret;
//}

//             T     -1 T   
// Compute   v1 *(J M  J )_mn * v2   where the indices mn are given by 
// the nodes associated with stations s1 & s2.
//
static double
computeA(const CDSVec3&    v1,
         const IVMLoopWNodes& loop1, int s1,
         const IVMLoopWNodes& loop2, int s2,
         const CDSVec3&    v2)
{
    const IVMRigidBodyNode* n1 = &loop1.tips(s1).getNode();
    const IVMRigidBodyNode* n2 = &loop2.tips(s2).getNode();

    CDSMat33 one(0.); one.setDiag(1.);

    CDSVec6 t1 = v1 * blockMat12(crossMat(n1->getOB_G() - loop1.tipPos(s1)),one);
    CDSVec6 t2 = blockMat21(crossMat(loop2.tipPos(s2) - n2->getOB_G()),one) * v2;

    while ( n1->getLevel() > n2->getLevel() ) {
        t1 = t1 * n1->getPsiT();
        n1 = n1->getParent();
    }

    while ( n2->getLevel() > n1->getLevel() ) {
        t2 = MatrixTools::transpose(n2->getPsiT()) * t2;
        n2 = n2->getParent();
    }

    while ( n1 != n2 ) {
        if (n1->isGroundNode() || n2->isGroundNode()  ) {
            t1.set(0.);  //not in same branch (or same tree -- sherm)
            t2.set(0.);
            cout << "computeA: cycles wasted calculating missed branch: "
                    << loop1.tips(s1) << " <-> " << loop2.tips(s2) << '\n';
            return 0.;
        }
        t1 = t1 * n1->getPsiT();
        t2 = MatrixTools::transpose(n2->getPsiT()) * t2;
        n1 = n1->getParent();
        n2 = n2->getParent();
    }

    // here n1==n2

    double ret = t1 * n1->getY() * t2;
    return ret;
}

//
// To be called for LengthSets consecutively from tip to base.
// This will calculate a force for every station in every loop
// contained in this IVMLengthSet, and store that force in the runtime
// block associated with that loop. It is up to the caller to do
// something with these forces.
//
// See Section 2.6 on p. 294 of Schwieters & Clore, 
// J. Magnetic Resonance 152:288-302. Equation reference below
// are to that paper. (sherm)
//
void
IVMLengthSet::calcConstraintForces() const
{
    // This is the acceleration error for each loop constraint in this
    // IVMLengthSet. We get a single scalar error per loop, since each
    // contains one distance constraint.
    // See Eq. [53] and the last term of Eq. [66].
    RVec0 rhs(loops.size(),0.);
    for (int i=0 ; i<loops.size() ; i++) {
        rhs(i) = abs2(loops[i].tipVel(2) - loops[i].tipVel(1))
                   + dot(loops[i].tipAcc(2) - loops[i].tipAcc(1) , 
                         loops[i].tipPos(2) - loops[i].tipPos(1));
    }

    // Here A = Q*(J inv(M) J')*Q' where J is the kinematic Jacobian for
    // the constrained points and Q is the constraint Jacobian. See first 
    // term of Eq. [66].
    RMat A(loops.size(),loops.size(),0.);
    for (int i=0 ; i<loops.size() ; i++) {
        const CDSVec3 v1 = loops[i].tipPos(2) - loops[i].tipPos(1);
        for (int bi=1 ; bi<=2 ; bi++)
            for (int bj=1 ; bj<=2 ; bj++) {
                double maxElem = 0.;
                for (int j=i ; j<loops.size() ; j++) {
                    const CDSVec3 v2 = loops[j].tipPos(2) - loops[j].tipPos(1);
                    double contrib = computeA(v1, loops[i], bi,
                                                  loops[j], bj, v2);
                    A(i,j) += contrib * (bi==bj ? 1 : -1);

                    if ( fabs(contrib) > maxElem ) maxElem = fabs(contrib);
                    if ( maxElem>0. && fabs(contrib)/maxElem < lConstraints->bandCut ) {
                        if ( verbose&IVMInternalDynamics::printLoopDebug )
                            cout << "IVMLengthSet::calcConstraintForces: setting A("
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

    for (int i=0 ; i<loops.size() ; i++)   //fill lower triangle
        for (int j=0 ; j<i ; j++)
            A(i,j) = A(j,i);

    //FIX: using inverse is inefficient
    const RVec0 lambda = inverse(A) * rhs;

    // add forces due to these constraints
    for (int i=0 ; i<loops.size() ; i++) {
        const CDSVec3 frc = lambda(i) * (loops[i].tipPos(2) - loops[i].tipPos(1));
        loops[i].setTipForce(2, -frc);
        loops[i].setTipForce(1,  frc);
    }
}

void IVMLengthSet::addInCorrectionForces(CDSVecVec6& spatialForces) const {
    for (int i=0; i<loops.size(); ++i) {
        for (int t=1; t<=2; ++t) {
            const IVMRigidBodyNode& node = loops[i].tipNode(t);
            const CDSVec3& force = loops[i].tipForce(t);
            const CDSVec3& moment = cross(loops[i].tipPos(t) - node.getOB_G(), force);
            spatialForces[node.getNodeNum()] += blockVec(moment, force);
        }
    }
}

void IVMLengthSet::testAccel()
{
    double testTol=1e-8;
    for (int i=0 ; i<loops.size() ; i++) {
        double test= dot( loops[i].tipAcc(2) - loops[i].tipAcc(1),
                          loops[i].tipPos(2) - loops[i].tipPos(1) ) +
                     abs2( loops[i].tipVel(2) - loops[i].tipVel(1) );
        if ( fabs(test) > testTol )
            cout << "IVMLengthSet::testAccel: constraint condition between atoms "
                 << loops[i].tips(1) << " and " << loops[i].tips(2) << " violated.\n"
                 << "\tnorm of violation: " << fabs(test) << '\n';
    }
    cout.flush();
}
   
//
// Fix internal force such that it obeys the constraint conditions.
//
void
IVMLengthSet::fixInternalForce(RVec& forceInternal)
{
    RMat  grad = calcGrad();
    if ( verbose & IVMInternalDynamics::printLoopDebug ) {
        RVec pos(rbTree.getDim()); rbTree.getPos(pos);
        testGrad(pos,grad);
    }
    RVec0 rhs  = multiForce(forceInternal,grad); 

    RMat A = MatrixTools::transpose(grad) * grad;

    //FIX: using inverse is inefficient
    RVec0 lambda = inverse(A) * rhs;

    // add forces due to these constraints
    addForce(forceInternal, grad * lambda); 
}

RVec0 
IVMLengthSet::multiForce(const RVec& forceInternal, const RMat& mat)
{
    RVec vec(dim);    //build vector same size as smallMat
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        int dim = nodeMap[j]->getDim();
        SubVec(vec,indx,dim) =
            ConstSubVec(forceInternal, nodeMap[j]->getStateOffset(),dim).vector();
        indx += dim;
    }

    //FIX: using transpose is inefficient
    RVec0 ret = MatrixTools::transpose(mat) * vec; 
    return ret;
}

void
IVMLengthSet::addForce(RVec& forceInternal,
                    const RVec& vec)
{
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        const int dim = nodeMap[j]->getDim();
        SubVec(forceInternal, nodeMap[j]->getStateOffset(),dim)
            -= ConstSubVec(vec,indx,dim);
        indx += dim;
    }
}

void
IVMLengthSet::testInternalForce(const RVec& forceInternal)
{
    RVec vec(dim);    //build vector same size as smallMat (index from 1)
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        const int dim = nodeMap[j]->getDim();
        SubVec(vec,indx,dim) = 
            ConstSubVec(forceInternal,nodeMap[j]->getStateOffset(),dim).vector();
        indx += dim;
    }

    RMat grad = calcGrad();
    RVec0 test = MatrixTools::transpose(grad) * vec;

    double testTol=1e-8;
    for (int i=0 ; i<loops.size() ; i++) {
        if ( test(i) > testTol )
            cout << "IVMLengthSet::Gradient: constraint condition between atoms "
                 << loops[i].tips(1) << " and " << loops[i].tips(2) << " violated.\n"
                 << "\tnorm of violation: " << test(i) << '\n';
    }
    cout.flush();
}

void
IVMLengthConstraints::fixVel0(RVec& iVel)
{
    // use molecule grouping
    for (int i=0 ; i<priv->accConstraints.size() ; i++)
        priv->accConstraints[i]->fixVel0(iVel);
}

//
// Correct internal velocities such that the constraint conditions are
// met under the condition that the disturbance to atom velocites is
// small as possible.
//
void
IVMLengthSet::fixVel0(RVec& iVel)
{
    // store internal velocities
    RVec iVel0 = iVel;

    CDSVector<double> rhs(loops.size());
    for (int k=0 ; k<loops.size() ; k++) 
        rhs(k) = dot(loops[k].tipPos(2) - loops[k].tipPos(1) ,
                     loops[k].tipVel(2) - loops[k].tipVel(1));

    RMat mat(loops.size(),loops.size());
    CDSList<RVec> deltaIVel(loops.size(),0);
    for (int m=0 ; m<loops.size() ; m++) {
        iVel.set(0.0);
        rbTree.setVel( iVel );
        deltaIVel[m].resize(iVel.size());

        // sherm: I think the following is a unit "probe" velocity, projected
        // along the separation vector. 
        // That would explain the fact that there are no velocities here!
        const CDSVec3 probeImpulse = loops[m].tipPos(2)-loops[m].tipPos(1);
        loops[m].setTipForce(2,  probeImpulse);
        loops[m].setTipForce(1, -probeImpulse);

        CDSVecVec6 spatialImpulse(rbTree.getNBodies());
        for (int ii=0; ii < rbTree.getNBodies(); ++ii)
            spatialImpulse[ii].set(0.);
        for (int t=1; t<=2; ++t) {
            const IVMRigidBodyNode& node = loops[m].tipNode(t);
            const CDSVec3& force = loops[m].tipForce(t);
            const CDSVec3& moment = cross(loops[m].tipPos(t) - node.getOB_G(), force);
            spatialImpulse[node.getNodeNum()] += blockVec(moment, force);
        }

        // calc deltaVa from k-th constraint condition
        //FIX: make sure that nodeMap is ordered by level
        // get all nodes in given molecule, ordered by level
        rbTree.calcZ(spatialImpulse);
        rbTree.calcTreeAccel();
        rbTree.getAcc(deltaIVel[m]);

        rbTree.setVel( deltaIVel[m] );

        for (int n=0 ; n<loops.size() ; n++)
            mat(n,m) = dot(loops[n].tipPos(2) - loops[n].tipPos(1),
                           loops[n].tipVel(2) - loops[n].tipVel(1));

        //store results of l-th constraint on deltaVa-k
    }
    RVec0 lambda = inverse(mat) * rhs;

    iVel = iVel0;
    for (l_int m=0 ; m<loops.size() ; m++)
        iVel -= lambda(m) * deltaIVel[m];

    rbTree.setVel(iVel);
}
   
// NOTHING BELOW HERE IS BEING USED
void 
IVMLengthSet::determineCouplings() 
{
 
// for (int i=0 ; i<loops.size() ; i++)
//   for (int j=0 ; j<loops.size() ; j++)
//     for (int bi=1 ; bi<=2 ; bi++)
//     for (int bj=1 ; bj<=2 ; bj++)
//       coupled[ loops[i].tips(bi)->node ][ loops[j].tips(bj)->node ] = 1;
 
// for (l_int i=0 ; i<loops.size() ; i++)
//   for (int j=0 ; j<loops.size() ; j++)
//     for (int bi=1 ; bi<=2 ; bi++)
//     for (IVMRigidBodyNode* node=loops[i].tips(bi)->node ;
//          node->levelA() ; 
//          node=node->parentA())
//       for (int bj=1 ; bj<=2 ; bj++)
//         if (node == loops[j].tips(bj)->node) 
//           coupled[ loops[i].tips(bi)->node ][ loops[j].tips(bj)->node ] = 1;
//
// // symmetrize
// for (Map_b2::iterator iti=coupled.begin() ; iti!=coupled.end() ; ++iti)
//   for (Map_b1::iterator itj=iti->second.begin() ; 
//      itj!=iti->second.end()                   ; 
//      ++itj                                    )
//     if ( itj->second )
//     coupled[ itj->first ][ iti->first ]  = 1;

       
// if ( IVMInternalDynamics::verbose&IVMInternalDynamics::printLoopDebug ) {
//   cout << "coupling matrix:\n";
//   for (l_int i=0 ; i<loops.size() ; i++)
//     for (int j=0 ; j<loops.size() ; j++)
//     for (int bi=1 ; bi<=2 ; bi++)
//       for (int bj=1 ; bj<=2 ; bj++)
//         cout << loops[i].tips(bi)->index << " " 
//          << loops[j].tips(bj)->index << ": "
//          << coupled[loops[i].tips(bi)->node][loops[j].tips(bj)->node] 
//          << '\n';
// }

}

//static bool
//sameBranch(const IVMRigidBodyNode* tip,
//           const IVMLoopWNodes& l )
//{
// for ( const IVMRigidBodyNode* node=tip ;
//     node->levelA()            ;
//     node=node->parentA()      )
//   if (node == l.tips(1)->node || node == l.tips(2)->node)
//     return 1;
// return 0;
//}
