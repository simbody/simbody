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

#include "dint-loop.h"
#include "dinternal.h"
#include "AtomTree.h"
#include "dint-node.h"
#include "dint-atom.h"
#include "vec3.h"
#include "cdsListAutoPtr.h"
#include "RVec.h"
#include <cdsMath.h>
#include <cdsMatrix.h>
#include <fixedVector.h>
#include <cdsVector.h>
#include <subVector.h>
#include <subMatrix.h>
#include <matrixTools.h>

#include <newtonRaphson.h>

#include <cdsIomanip.h>

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using MatrixTools::inverse;
#endif /* USE_CDS_NAMESPACE */

typedef CDSMatrix<double>     RMat;
typedef CDSVector<double,0>   RVec0;
typedef SubVector<RVec>       SubVec;
typedef SubVector<const RVec> ConstSubVec;
typedef SubMatrix<RMat>       SubMat;

static int compareLevel(const LoopWNodes& l1,    //forward declarations
                        const LoopWNodes& l2);

//static bool sameBranch(const HingeNode* tip,
//                       const LoopWNodes& l );

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
class LoopWNodes {
    HingeNode*                base; // highest-level common ancestor of tips
    FixedVector<IVMAtom*,2,1> tips; // the two connected atoms, sorted so that
                                    //   tips(1).level <= tips(2).level
    FixedVector<NodeList,2,1> nodes;// the two paths: base..tip1, base..tip2,
                                    //   incl. tip nodes but not base
    const HingeNode*          moleculeNode;
public:
    LoopWNodes() {}
    LoopWNodes(const Loop& l);

    friend class LengthSet;
    friend ostream& operator<<(ostream& os, const LengthSet& s);
    friend void LengthConstraints::construct(CDSList<Loop>&);
    friend int compareLevel(const LoopWNodes& l1,
                            const LoopWNodes& l2);
    friend bool sameBranch(const HingeNode* tip,
                           const LoopWNodes& l );
};

LoopWNodes::LoopWNodes(const Loop& l) 
{
    using InternalDynamics::Exception;


    if ( l.tip1->node == l.tip2->node ) {
        cout << "LoopWNodes::LoopWNodes: bad topology:\n\t"
             << "loop atoms " << l.tip1
             << " and  "      << l.tip2
             << " are now in the same node. Deleting loop.\n";
        throw BadNodeDef();
    }

    // Ensure that tips(2) is the atom which is farther from the base.
    if ( l.tip1->node->getLevel() > l.tip2->node->getLevel() )
         tips(1) = l.tip2, tips(2) = l.tip1;
    else tips(1) = l.tip1, tips(2) = l.tip2;

    // Collect up the node path from tips(2) down to the last node on its
    // side of the loop which is at a higher level than tips(1) (may be none).
    HingeNode* node1 = tips(1)->node;
    HingeNode* node2 = tips(2)->node;
    while ( node2->getLevel() > node1->getLevel() ) {
        nodes(2).append(node2);
        node2 = node2->getParent();
    }

    // We're at the same level on both sides of the loop. Run down both
    // simultaneously until we hit the first common ancestor, collecting
    // up the nodes along the two paths (but not the common ancestor).
    while ( node1 != node2 ) {
        if ( node1->isGroundNode() ) {
            cerr << "LoopWNodes::LoopWNodes: could not find base node.\n\t"
                 << "loop between atoms " << tips(1)->index << " and " 
                 << tips(2)->index << "\n";
            throw Exception("LoopWNodes::LoopWNodes: could not find base atom");
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
        cerr << "LoopWNodes::LoopWNodes: could not find molecule node.\n\t"
             << "loop between atoms " << tips(1)->index << " and " 
             << tips(2)->index << "\n";
        throw Exception("LoopWNodes::LoopWNodes: could not find molecule node");
    }
}

typedef CDSList<LoopWNodes> LoopList;

class LengthSet {
    static void construct(const LoopList& loops);
    const LengthConstraints* lConstraints;
    IVM*                     ivm;
    LoopList                 loops;    
    CDSList<double>          lengths;
    int                      dim;
    CDSList<HingeNode*>      nodeMap; //unique nodes (intersection of loops->nodes)
public:
    LengthSet(const LengthConstraints* lConstraints, const LoopWNodes& loop)
        : lConstraints(lConstraints), ivm(lConstraints->ivm), dim(0)
    {
        addConstraint(loop);
    }

    //  ~LengthSet() { for (int i=0 ; i<loops.size() ; i++) delete loops[i]; }

    void addConstraint(const LoopWNodes& loop) {
        loops.append( LoopWNodes(loop) );
        lengths.append( sqrt(abs2(loop.tips(1)->pos - loop.tips(2)->pos)) );
        for (int b=1 ; b<=2 ; b++)
            for (int i=0 ; i<loop.nodes(b).size() ; i++)
                if ( !nodeMap.contains(loop.nodes(b)[i]) ) {
                    dim += loop.nodes(b)[i]->getDim();
                    nodeMap.append( loop.nodes(b)[i] );
                }
    }

    void determineCouplings();
    bool contains(HingeNode* node) {
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

    void  setPosVel(const RVec& pos, const RVec& vel) const;
    RVec  getPos();
    RVec0 calcPosB(const RVec&, const RVec&) const;
    RVec0 calcVelB(const RVec&, const RVec&) const;
    RVec  calcPosZ(const RVec& b) const;
    RMat  calcGrad() const;
    RMat  calcGInverse() const;
    void  fixAccel();
    void  fixVel0(RVec&);
    void  fixInternalForce(RVec&);
    RVec0 multiForce(const RVec&, const RMat& mat);
    void  addForce(RVec&, const RVec& ve);


    //  void  enforce(RVec& pos, RVec& vel);
    void testAccel();
    void testInternalForce(const RVec&);
    friend ostream& operator<<(ostream& os, const LengthSet& s);
    friend class CalcPosB;
    friend class CalcVelB;
    friend class CalcPosZ;
    friend class CalcVelZ;

    void fdgradf(const RVec& pos, const RVec& vel, RMat& grad) const;
    void testGrad(const RVec& pos, const RMat& grad) const;
};

ostream& 
operator<<(ostream& os, const LengthSet& s) 
{
    for (int i=0 ; i<s.loops.size() ; i++)
        os << setw(4) << s.loops[i].tips(1) << "<->" 
           << setw(4) << s.loops[i].tips(2)  << ": " 
           << s.lengths[i] << "  ";
    return os;
}


class LengthConstraintsPrivates {
public:
    CDSListAutoPtr<LengthSet> constraints;     // used for pos, vel
    CDSListAutoPtr<LengthSet> accConstraints;  // used for acc
    NewtonRaphson posMin, velMin;
    LengthConstraintsPrivates() : posMin(cout), velMin(cout) {}
};


LengthConstraints::LengthConstraints(IVM* ivm)
    : bandCut(1e-7), maxIters( 20 ), maxMin( 20 ), ivm(ivm)
{
    priv = new LengthConstraintsPrivates;
    priv->posMin.maxIters = maxIters;
    priv->posMin.maxMin   = maxMin;
    priv->posMin.tol      = ivm->Ctolerance();
    priv->velMin.maxIters = maxIters;
    priv->velMin.maxMin   = maxMin;
    priv->velMin.tol      = ivm->Ctolerance()*10;
}

LengthConstraints::~LengthConstraints()
{
    delete priv;
}

//template<class CalcB,class CalcZ>
//void
//NewtonRaphson(RVec& x,
//        CalcB calcB,
//        CalcZ calcZ,
//        const double& tol);

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

void 
LengthSet::determineCouplings() 
{
 
// for (int i=0 ; i<loops.size() ; i++)
//   for (int j=0 ; j<loops.size() ; j++)
//     for (int bi=1 ; bi<=2 ; bi++)
//     for (int bj=1 ; bj<=2 ; bj++)
//       coupled[ loops[i].tips(bi)->node ][ loops[j].tips(bj)->node ] = 1;
 
// for (l_int i=0 ; i<loops.size() ; i++)
//   for (int j=0 ; j<loops.size() ; j++)
//     for (int bi=1 ; bi<=2 ; bi++)
//     for (HingeNode* node=loops[i].tips(bi)->node ;
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

       
// if ( InternalDynamics::verbose&InternalDynamics::printLoopDebug ) {
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
//sameBranch(const HingeNode* tip,
//         const LoopWNodes& l )
//{
// for ( const HingeNode* node=tip ;
//     node->levelA()            ;
//     node=node->parentA()      )
//   if (node == l.tips(1)->node || node == l.tips(2)->node)
//     return 1;
// return 0;
//}

//1) construct: given list of loops 
//   a) if appropriate issue warning and exit.
//   b) sort by node 1 of each loop.
//   c) find loops which intersect: combine loops and increment
//    number of length constraints
void
LengthConstraints::construct(CDSList<Loop>& iloops)
{
    //clean up
    priv->constraints.resize(0);
    priv->accConstraints.resize(0);

    // if ( !useLengthConstraint )  FIX:!!!
    //   //issue error message?
    //   return;

    LoopList loops;
    for (int i=0 ; i<iloops.size() ; i++) {
        try {
            LoopWNodes loop( iloops[i] );
            loops.append( loop );
        }
        catch ( BadNodeDef ) {}
    }

    // sort loops by nodes[0]->level
    loops.sort(compareLevel);
    LoopList accLoops = loops;  //version for acceleration

    // find intersections -- this version keeps hierarchical loops distinct
    for (l_int i=0 ; i<loops.size() ; i++) {
        priv->constraints.append( new LengthSet(this, loops[i]) );
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
    for (l_int i=0 ; i<accLoops.size() ; i++) {
        priv->accConstraints.append( new LengthSet(this, accLoops[i]) );
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

    // int compareBaseNodeLevel(const LengthSet& l1,
    //                const LengthSet& l2) {
    //   return ( (l1->baseAtom() > l2->baseAtom())?1 :
    //          ( (l1->baseAtom()==l2->baseAtom())?0 : -1) );
    // constraints.sort(compareBaseNodeLevel);

    if (priv->constraints.size()>0 && ivm->verbose()&InternalDynamics::printLoopInfo) {
        cout << "LengthConstraints::construct: pos/vel length constraints found:\n";
        for (l_int i=0 ; i<priv->constraints.size() ; i++)
            cout << "\t" << *priv->constraints[i] << "\n";
        cout << "LengthConstraints::construct: accel length constraints found:\n";
        for (l_int i=0 ; i<priv->accConstraints.size() ; i++)
            cout << "\t" << *priv->accConstraints[i] << "\n";
    }
}

class CalcPosB {
    const LengthSet* constraint;
    const RVec& vel;
public:
    CalcPosB(const LengthSet* constraint, const RVec& vel)
        : constraint(constraint), vel(vel) {}
    RVec0 operator()(const RVec& pos) 
        { return constraint->calcPosB(pos,vel); }
};

class CalcPosZ {
    const LengthSet* constraint;
public:
    CalcPosZ(const LengthSet* constraint) : constraint(constraint) {}
    RVec0 operator()(const RVec& b) 
        { return constraint->calcPosZ(b); }
};

class CalcVelB {
    const LengthSet* constraint;
    const RVec& pos;
public:
    CalcVelB(const LengthSet* constraint, const RVec& pos)
        : constraint(constraint), pos(pos) {}
    RVec0 operator()(const RVec& vel) 
        { return constraint->calcVelB(pos,vel); }
};

//
// calculate the constraint violation (zero when constraint met)
//
RVec0
LengthSet::calcPosB(const RVec& pos, const RVec& vel) const
{
    setPosVel(pos,vel);

    RVec0 b( loops.size() );
    for (int i=0 ; i<loops.size() ; i++) 
        b(i) = lengths[i] - 
                sqrt(abs2(loops[i].tips(1)->pos - loops[i].tips(2)->pos));
    return b;
}

//
// calculate the constraint violation (zero when constraint met)
//
RVec0
LengthSet::calcVelB(const RVec& pos, const RVec& vel) const 
{
    setPosVel(pos,vel);

    RVec0 b( loops.size() );
    for (int i=0 ; i<loops.size() ; i++) {
        const IVMAtom* tip1 = loops[i].tips(1);
        const IVMAtom* tip2 = loops[i].tips(2);
        b(i) = -dot( unitVec(tip2->pos - tip1->pos) , tip2->vel - tip1->vel );
    }

    return b;
}

//
// Given a vector containing violations of position constraints, calculate
// a state update which should drive those violations to zero (if they
// were linear).
//
RVec
LengthSet::calcPosZ(const RVec& b) const
{
    RVec dir = calcGInverse() * b;
    RVec z(ivm->dim(),0.0);

    // map the vector dir back to the appropriate elements of z
    for (int i=0 ; i<nodeMap.size() ; i++) 
        SubVec(z,nodeMap[i]->offset(), nodeMap[i]->getDim())
            = SubVec(dir,i+1,nodeMap[i]->getDim());

    return z;
}

//
// Given a vector containing violations of velocity constraints, calculate
// a state update which would drive those violations to zero (if they
// are linear and well conditioned).
//
class CalcVelZ {
    const LengthSet* constraint;
    const RMat       gInverse;
public:
    CalcVelZ(const LengthSet* constraint) : 
    constraint(constraint), gInverse(constraint->calcGInverse()) {}
    RVec operator()(const RVec& b) {
        RVec dir = gInverse * b;
        RVec z(constraint->ivm->dim(),0.0);

        // map the vector dir back to the appropriate elements of z
        for (int i=0 ; i<constraint->nodeMap.size() ; i++) {
            const HingeNode* n = constraint->nodeMap[i];
            SubVec(z,n->offset(),n->getDim()) =  SubVec(dir,i+1,n->getDim());
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
LengthConstraints::enforce(RVec& pos, RVec& vel)
{
    priv->posMin.verbose  = ((ivm->verbose()&InternalDynamics::printLoopDebug) != 0);
    priv->velMin.verbose  = ((ivm->verbose()&InternalDynamics::printLoopDebug) != 0);
    try { 
        for (int i=0 ; i<priv->constraints.size() ; i++) {
            if ( ivm->verbose()&InternalDynamics::printLoopDebug )
                cout << "LengthConstraints::enforce: position " 
                     << *(priv->constraints[i]) << '\n';
            priv->posMin.calc(pos,
                              CalcPosB(priv->constraints[i].get(),vel),
                              CalcPosZ(priv->constraints[i].get()));
        }
        for (l_int i=0 ; i<priv->constraints.size() ; i++) {
            if ( ivm->verbose()&InternalDynamics::printLoopDebug )
                cout << "LengthConstraints::enforce: velocity " 
                     << *priv->constraints[i] << '\n';
            priv->velMin.calc(vel,
            CalcVelB(priv->constraints[i].get(),pos),
            CalcVelZ(priv->constraints[i].get()));
        }
    }
    catch ( NewtonRaphson::Fail cptn ) {
        cout << "LengthConstraints::enforce: exception: "
             << cptn.mess << '\n';
    } 
// catch ( ... ) {
//   cout << "LengthConstraints::enforce: exception: "
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
//
//
//
void
LengthSet::setPosVel(const RVec& pos,
             const RVec& vel) const
{
    for (int i=0 ; i<nodeMap.size() ; i++)
        nodeMap[i]->setPosVel(pos,vel); //also calc necessary properties
}

//
////A = df / dtheta_i
////  = [ -(q_ni-qn0)x , 1 ] [ phi_na^T Ha , phi_nb^T Hb , ... ]^T
////
//

// Presumes that calcEnergy has been called previously with current 
// value of ipos.
void 
LengthSet::fdgradf( const RVec& pos,
                    const RVec& vel,
                    RMat&       grad) const 
{
    // Gradf gradf(tree);
    // gradf(x,grad); return;
    double eps = 1e-8;
    CalcPosB calcB(this,vel);
    RVec0 b = calcB(pos);
    int grad_indx=0;
    for (int i=0 ; i<nodeMap.size() ; i++) {
        int pos_indx=nodeMap[i]->offset();
        for (int j=0 ; j<nodeMap[i]->getDim() ; j++,pos_indx++,grad_indx++) {
            RVec posp = pos;
            posp(pos_indx) += eps;
            RVec0 bp = calcB(posp);
            for (int k=0 ; k<b.size() ; k++)
                grad(grad_indx,k) = -(bp(k)-b(k)) / eps;
        }
    }
}

void 
LengthSet::testGrad(const RVec& pos, const RMat& grad) const
{
    double tol = 1e-4;
    // Costf costf(tree);
    // double f = costf(x);
    RMat fdgrad(dim,loops.size());
    RVec vel(pos.size(),0.);
    fdgradf(pos,vel,fdgrad);

    for (int i=0 ; i<grad.rows() ; i++)
        for (int j=0 ; j<grad.cols() ; j++)
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
RMat
LengthSet::calcGrad() const
{
    RMat grad(dim,loops.size(),0.0);
    Mat33 one(0.0); one.setDiag(1.0);  //FIX: should be done once
    for (int i=0 ; i<loops.size() ; i++) {
        const LoopWNodes& l = loops[i];
        FixedVector<CDSList<Mat66>,2,1> phiT;
        for (int b=1 ; b<=2 ; b++) {
            phiT(b).resize( l.nodes(b).size() );
            if ( l.nodes(b).size() ) {
                phiT(b)[phiT(b).size()-1].set(0.0);
                phiT(b)[phiT(b).size()-1].setDiag(1.0);
                for (int j=l.nodes(b).size()-2 ; j>=0 ; j-- ) {
                    HingeNode* n = l.nodes(b)[j+1];
                    phiT(b)[j] = phiT(b)[j+1] * transpose( n->getPhi() );
                }
            }
        }
        // compute gradient
        Vec3 uBond = unitVec(l.tips(2)->pos - l.tips(1)->pos);
        FixedVector<FixedMatrix<double,3,6>,2,1> J;
        for (l_int b=1 ; b<=2 ; b++)
            J(b) = blockMat12(-crossMat(l.tips(b)->pos -
                              l.tips(b)->node->getAtom(0)->pos),one);
        int g_indx=0;
        for (int j=0 ; j<nodeMap.size() ; j++) {
            RMat H = MatrixTools::transpose(nodeMap[j]->getH());
            double elem=0.0;
            int l1_indx = l.nodes(1).getIndex(nodeMap[j]);
            int l2_indx = l.nodes(2).getIndex(nodeMap[j]);
            for (int k=0 ; k<H.cols() ; k++) {
                Vec6 Hcol = subCol(H,k,0,5).vector();
                if ( l1_indx >= 0 ) { 
                    elem = -dot(uBond , Vec3(J(1) * phiT(1)[l1_indx]*Hcol));
                } else if ( l2_indx >= 0 ) { 
                    elem = dot(uBond , Vec3(J(2) * phiT(2)[l2_indx]*Hcol));
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
    for (int i=0 ; i<m.rows() ; i++)
        for (int j=0 ; j<m.cols() ; j++)
            ret += CDSMath::sq( m(i,j) );
    return ret;
}

//
// calculate generalized inverse which minimizes changes in soln vector
//
RMat
LengthSet::calcGInverse() const
{
    RMat grad = calcGrad();
    if ( ivm->verbose() & InternalDynamics::printLoopDebug )
        testGrad(ivm->tree()->getPos(),grad);

    RMat ret(grad.rows(),grad.rows(),0.0);
    if ( abs2(grad) > 1e-10 ) 
        ret = grad * inverse( MatrixTools::transpose(grad)*grad );
    return ret;
}

//template<class CalcB,class CalcZ>
//void
//NewtonRaphson(RVec& x,
//        CalcB calcB,
//        CalcZ calcZ,
//        const double& tol) 
//{
// using InternalDynamics::Exception;
//
// RVec0 b = calcB(x);
// double norm = sqrt(abs2(b)) / x.size();
// double onorm=norm;
// int  iters=0;
// cout.setf(ios::fixed);
// if ( ivm->verbose&InternalDynamics::printLoopDebug )
//   cout << "LengthSet::enforce: iter: " 
//      << iters << "  norm: " << norm << '\n';
// cout.unsetf(ios::fixed);
// bool finished=(norm < tol);
// while (!finished) {
//   RVec ox = x;
//   
//   RVec z = calcZ(b);
//   //   z *= gradStep; //gradStep should be 1
//   x += z;
//   
//   iters++;
//   
//   b = calcB(x);
//   norm = sqrt(abs2(b)) / x.size();
//   
//   if ( norm > onorm &&
//      verbose&InternalDynamics::printLoopDebug)
//     cout << "LengthSet::enforce: newton-Raphson failed."
//        << " Trying gradient search.\n";
//   
//   int mincnt = maxMin;
//   z *= -1.0;
//   while ( norm > onorm ) {
//     if ( mincnt < 1 ) 
//       throw Exception("LengthSet::enforce: too many minimization steps taken.\n");
//     z *= 0.5;
//     x = ox + z;
//     b = calcB(x);
//     norm = sqrt(abs2(b)) / x.size();
//     if ( verbose&InternalDynamics::printLoopDebug )
//     cout << "LengthSet::enforce: iter: " 
//          << iters << "  norm: " << norm << '\n';
//     mincnt--;
//   }
//
//   onorm = norm;
//   
//   if ( verbose&InternalDynamics::printLoopDebug )
//     cout << "LengthSet::enforce: iter: " 
//        << iters << "  norm: " << norm << '\n';
//     
//   if (norm < tol)
//     finished=1;
//   if (iters > maxIters) 
//     throw Exception("LengthSet::enforce: maxiters exceeded");
// }
//} /* LengthConstraints::NewtonRaphson */

//acceleration:
//  0) after initial acceleration calculation:
//  1) calculate Y (block diagonal matrix) for all nodes (base to tip)
//  2) for each LengthSet, calculate Lagrange multiplier(s) and add in
//     resulting force- update accelerations
//     Do this step from tip to base.
//

bool 
LengthConstraints::fixAccel()
{
    if ( priv->accConstraints.size() == 0 )
        return 0;

    ivm->calcY();

    for (int i=priv->accConstraints.size()-1 ; i>=0 ; i--)
        priv->accConstraints[i]->fixAccel();

    if ( ivm->verbose()&InternalDynamics::printLoopInfo )
        for (l_int i=priv->accConstraints.size()-1 ; i>=0 ; i--) 
            priv->accConstraints[i]->testAccel();

    return 1;
}

void 
LengthConstraints::fixGradient(RVec& forceInternal)
{
    if ( priv->constraints.size() == 0 )
        return;

    for (int i=priv->constraints.size()-1 ; i>=0 ; i--)
        priv->constraints[i]->fixInternalForce(forceInternal);

    for (l_int i=priv->constraints.size()-1 ; i>=0 ; i--) 
        priv->constraints[i]->testInternalForce(forceInternal);
}

typedef SubVector<const Vec6>       SubVec6;
typedef SubVector<const Vec6>       SubVec6;

//
// Calculate the acceleration of atom assuming that the spatial acceleration
// of the body (node) it's on is available.
//
static Vec3
getAccel(const IVMAtom* a)
{
    const HingeNode* n = a->node;
    Vec3 ret( SubVec6(n->getSpatialAcc(),3,3).vector() );
    ret += cross( SubVec6(n->getSpatialAcc(),0,3).vector() , a->pos - n->getAtom(0)->pos );
    ret += cross( SubVec6(n->getSpatialVel(),0,3).vector() , a->vel - n->getAtom(0)->vel );
    return ret;
}

//             T     -1 T   
// Compute   v1 *(J M  J )_mn * v2   where the indices mn are given by 
// the nodes associated with atoms a1 & a2.
//
static double
computeA(const Vec3&    v1,
         const IVMAtom* a1,
         const IVMAtom* a2,
         const Vec3&    v2)
{
    const HingeNode* n1 = a1->node;
    const HingeNode* n2 = a2->node;

    Mat33 one(0.); one.setDiag(1.);

    Vec6 t1 = v1 * blockMat12(crossMat(n1->getAtom(0)->pos - a1->pos),one);
    Vec6 t2 = blockMat21(crossMat(a2->pos - n2->getAtom(0)->pos),one) * v2;

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
                    << a1 << " <-> " << a2 << '\n';
            return 0.;
        }
        t1 = t1 * n1->getPsiT();
        t2 = MatrixTools::transpose(n2->getPsiT()) * t2;
        n1 = n1->getParent();
        n2 = n2->getParent();
    }

    // here n1==n2

    return t1 * n1->getY() * t2;
}

//
// To be called for LengthSets consecutively from tip to base.
//
// See Section 2.6 on p. 294 of Schwieters & Clore, 
// J. Magnetic Resonance 152:288-302. Equation reference below
// are to that paper. (sherm)
//
void
LengthSet::fixAccel()
{
    // This is the acceleration error for each loop constraint in this
    // LengthSet. We get a single scalar error per loop, since each
    // contains one distance constraint.
    // See Eq. [53] and the last term of Eq. [66].
    RVec0 rhs(loops.size(),0.);
    for (int i=0 ; i<loops.size() ; i++) {
        const IVMAtom* t = loops[i].tips(2);
        const IVMAtom* b = loops[i].tips(1);
        rhs(i) = -(abs2(t->vel - b->vel)
                   + dot(getAccel(t) - getAccel(b) , t->pos - b->pos));
    }

    // Here A = Q*(J inv(M) J')*Q' where J is the kinematic Jacobian for
    // the constrained points and Q is the constraint Jacobian. See first 
    // term of Eq. [66].
    RMat A(loops.size(),loops.size(),0.);
    for (l_int i=0 ; i<loops.size() ; i++) {
        const Vec3 v1 = loops[i].tips(2)->pos - loops[i].tips(1)->pos;
        for (int bi=1 ; bi<=2 ; bi++)
            for (int bj=1 ; bj<=2 ; bj++) {
                double maxElem = 0.;
                for (int j=i ; j<loops.size() ; j++) {
                    const Vec3 v2 = loops[j].tips(2)->pos - loops[j].tips(1)->pos;
                    double contrib = computeA(v1, loops[i].tips(bi),
                                              loops[j].tips(bj), v2);
                    A(i,j) += contrib * (bi==bj ? 1 : -1);

                    if ( fabs(contrib) > maxElem ) maxElem = fabs(contrib);
                    if ( maxElem>0. && fabs(contrib)/maxElem < lConstraints->bandCut ) {
                        if ( ivm->verbose()&InternalDynamics::printLoopDebug )
                            cout << "LengthSet::fixAccel: setting A("
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

    for (l_int i=0 ; i<loops.size() ; i++)   //fill lower triangle
        for (int j=0 ; j<i ; j++)
            A(i,j) = A(j,i);

    //FIX: using inverse is inefficient
    const RVec0 lambda = inverse(A) * rhs;

    // add forces due to these constraints
    for (l_int i=0 ; i<loops.size() ; i++) {
        IVMAtom* b = loops[i].tips(1);
        IVMAtom* t = loops[i].tips(2);
        const Vec3 frc = lambda(i) * (t->pos - b->pos);
        t->deriv -= frc;
        b->deriv += frc;
    }

    ivm->updateAccel();
}

void
LengthSet::testAccel()
{
    double testTol=1e-8;
    for (int i=0 ; i<loops.size() ; i++) {
        IVMAtom* t = loops[i].tips(2);
        IVMAtom* b = loops[i].tips(1);
        double test= dot( getAccel(t) - getAccel(b) , t->pos - b->pos ) +
                     abs2( t->vel - b->vel );
        if ( fabs(test) > testTol )
            cout << "LengthSet::testAccel: constraint condition between atoms "
                 << b << " and " << t << " violated.\n"
                 << "\tnorm of violation: " << fabs(test) << '\n';
    }
    cout.flush();
}
   
//
// Fix internal force such that it obeys the constraint conditions.
//
void
LengthSet::fixInternalForce(RVec& forceInternal)
{
    RMat  grad = calcGrad();
    if ( ivm->verbose() & InternalDynamics::printLoopDebug )
        testGrad(ivm->tree()->getPos(),grad);
    RVec0 rhs  = multiForce(forceInternal,grad); 

    RMat A = MatrixTools::transpose(grad) * grad;

    //FIX: using inverse is inefficient
    RVec0 lambda = inverse(A) * rhs;

    // add forces due to these constraints
    addForce(forceInternal, grad * lambda); 
}

RVec0 
LengthSet::multiForce(const RVec& forceInternal, const RMat& mat)
{
    RVec vec(dim);    //build vector same size as smallMat
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        int dim = nodeMap[j]->getDim();
        SubVec(vec,indx,dim) =
            ConstSubVec(forceInternal, nodeMap[j]->offset(),dim).vector();
        indx += dim;
    }

    //FIX: using transpose is inefficient
    RVec0 ret = MatrixTools::transpose(mat) * vec; 
    return ret;
}

void
LengthSet::addForce(      RVec& forceInternal,
            const RVec& vec)
{
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        int dim = nodeMap[j]->getDim();

        SubVec(forceInternal, nodeMap[j]->offset(),dim)
            -= ConstSubVec(vec,indx,dim);
        indx += nodeMap[j]->getDim();
    }
}

void
LengthSet::testInternalForce(const RVec& forceInternal)
{
    RVec vec(dim);    //build vector same size as smallMat
    int indx=1;
    for (int j=0 ; j<nodeMap.size() ; j++) {
        int dim = nodeMap[j]->getDim();
        SubVec(vec,indx,dim) = 
        ConstSubVec(forceInternal,nodeMap[j]->offset(),dim).vector();
        indx += nodeMap[j]->getDim();
    }

    RMat grad = calcGrad();
    RVec0 test = MatrixTools::transpose(grad) * vec;

    double testTol=1e-8;
    for (int i=0 ; i<loops.size() ; i++) {
        IVMAtom* t = loops[i].tips(2);
        IVMAtom* b = loops[i].tips(1);
        if ( test(i) > testTol )
            cout << "LengthSet::Gradient: constraint condition between atoms "
                 << b << " and " << t << " violated.\n"
                 << "\tnorm of violation: " << test(i) << '\n';
    }
    cout.flush();
}

void
LengthConstraints::fixVel0(RVec& iVel)
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
LengthSet::fixVel0(RVec& iVel)
{
    // store internal velocities
    RVec iVel0 = iVel;

    CDSVector<double> rhs(loops.size());
    for (int k=0 ; k<loops.size() ; k++) 
        rhs(k) = dot(loops[k].tips(2)->pos - loops[k].tips(1)->pos ,
                     loops[k].tips(2)->vel - loops[k].tips(1)->vel);

    RMat mat(loops.size(),loops.size());
    CDSList<RVec> deltaIVel(loops.size(),0);
    for (int m=0 ; m<loops.size() ; m++) {
        iVel.set(0.0);
        ivm->tree()->setVel( iVel );
        loops[m].tips(2)->vel  = (loops[m].tips(2)->pos-loops[m].tips(1)->pos);
        loops[m].tips(1)->vel  = -loops[m].tips(2)->vel / loops[m].tips(1)->mass;
        loops[m].tips(2)->vel /=  loops[m].tips(2)->mass;

        // calc deltaVa from k-th constraint condition
        //FIX: make sure that nodeMap is ordered by level
        // get all nodes in given molecule, ordered by level
        ivm->tree()->propagateSVel();
        ivm->tree()->calcZ();

        // make sure that internal velocities of each node are correct

        deltaIVel[m] = ivm->tree()->calcGetAccel();
        ivm->tree()->setVel( deltaIVel[m] );

        for (int n=0 ; n<loops.size() ; n++)
            mat(n,m) = dot(loops[n].tips(2)->pos - loops[n].tips(1)->pos,
                           loops[n].tips(2)->vel - loops[n].tips(1)->vel);

        //store results of l-th constraint on deltaVa-k
    }
    RVec0 lambda = inverse(mat) * rhs;

    iVel = iVel0;
    for (l_int m=0 ; m<loops.size() ; m++)
        iVel -= lambda(m) * deltaIVel[m];

    ivm->tree()->setVel(iVel);
}
   
