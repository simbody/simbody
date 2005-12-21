/**@file
 *
 * Implementation for atom- and bond-sensitive model building.
 */

#include "AtomTree.h"
#include "internalDynamics.h"
#include "dinternal.h"

#include "AtomClusterNode.h"
#include "RigidBodyNode.h"

#include "dint-atom.h"
#include "cdsVec3.h"

#include <sthead.h>
#include <cdsMath.h>
#include <cdsString.h>
#include <cdsSStream.h>
#include <cdsVector.h>
#include <fixedVector.h>
#include <subVector.h>
#include <fixedMatrix.h>
#include <matrixTools.h>
#include <cdsIomanip.h>
#include <cdsFstream.h>

using namespace InternalDynamics;
using MatrixTools::inverse;

typedef FixedVector<double,6> CDSVec6;
typedef FixedMatrix<double,6> CDSMat66;


//
// structure of atomtree:
//  doubly linked tree
//

/**
 * Temporary class used in building molecule trees.
 */
class AT_Build {
public:
    IVM* ivm;
    CDSList<AtomLoop>& loops;
    AtomTree* const tree;
    CDSVector<bool,0> assignedAtoms; //used in locating cycle links
    AT_Build(IVM*,
             AtomTree* const, 
             AtomClusterNode*,
             CDSList<AtomLoop>&);
    void buildAtomClusterNode(AtomClusterNode*);
};

AtomTree::~AtomTree() {
    for (int i=0; i<nodeTree.size(); i++) {
        for (int j=0 ; j<nodeTree[i].size() ; j++) 
            delete nodeTree[i][j];
        nodeTree[i].resize(0);
    }
    nodeTree.resize(0);
}

void AtomTree::setClusterVelFromSVel(int level, int indx, const CDSVec6& sVel) {
    AtomClusterNode& ac = *nodeTree[level][indx];
    RigidBodyNode& rb = rbTree.updRigidBodyNode(ac.getRBIndex());
    return rb.setVelFromSVel(sVel);
}

void AtomTree::calcAtomPos() {
    for (int l=0; l<nodeTree.size(); l++)
        for (int j=0; j<nodeTree[l].size(); j++) {
            AtomClusterNode&     ac = *nodeTree[l][j];
            const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
            ac.calcAtomPos(rb.getR_GB(), rb.getOB_G());
        }
}

void AtomTree::calcAtomVel() {
    for (int l=0; l<nodeTree.size(); l++)
        for (int j=0; j<nodeTree[l].size(); j++) {
            AtomClusterNode&     ac = *nodeTree[l][j];
            const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
            ac.calcAtomVel(rb.getSpatialVel());
        }
}

void AtomTree::calcSpatialForces() {
    spatialForces.resize(rbTree.getNBodies());
    for (int l=0; l<nodeTree.size(); l++)
        for (int j=0; j<nodeTree[l].size(); j++) {
            AtomClusterNode&     ac = *nodeTree[l][j];
            ac.calcSpatialForce(spatialForces[ac.getRBIndex()]);
        }
}

// note that results go into spatialForces vector
void AtomTree::calcSpatialImpulses() {
    spatialForces.resize(rbTree.getNBodies());
    for (int l=0; l<nodeTree.size(); l++)
        for (int j=0; j<nodeTree[l].size(); j++) {
            AtomClusterNode&     ac = *nodeTree[l][j];
            ac.calcSpatialImpulse(spatialForces[ac.getRBIndex()]);
        }
}


double AtomTree::getClusterMass(int level, int indx) const {
    const AtomClusterNode& ac = *nodeTree[level][indx];
    const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
    return rb.getMass();
}

const CDSVec3& AtomTree::getClusterCOM_G(int level, int indx) const {
    const AtomClusterNode& ac = *nodeTree[level][indx];
    const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
    return rb.getCOM_G();
}

const CDSVec6& AtomTree::getClusterSpatialVel(int level, int indx) const {
    const AtomClusterNode& ac = *nodeTree[level][indx];
    const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
    return rb.getSpatialVel();
}

double AtomTree::calcClusterKineticEnergy(int level, int indx) const {
    const AtomClusterNode& ac = *nodeTree[level][indx];
    const RigidBodyNode& rb = rbTree.getRigidBodyNode(ac.getRBIndex());
    return rb.calcKineticEnergy();
}

RVec AtomTree::calcGetAccel() { 
    RVec acc(getIVMDim()); 
    rbTree.calcTreeAccel();
    rbTree.getAcc(acc);
    return acc;
}

void AtomTree::fixVel0(RVec& vel) {
    rbTree.fixVel0(vel);
}

// One-stop shopping -- calculates everything except the
// position & velocity kinematics.
RVec AtomTree::getAccel() {
    calcSpatialForces();
    rbTree.prepareForDynamics();
    rbTree.calcLoopForwardDynamics(spatialForces);
    RVec acc(getIVMDim()); 
    rbTree.getAcc(acc);
    return acc;
}

RVec AtomTree::getAccelIgnoringConstraints() {
    calcSpatialForces();
    rbTree.prepareForDynamics();
    rbTree.calcTreeForwardDynamics(spatialForces);
    RVec acc(getIVMDim()); 
    rbTree.getAcc(acc);
    return acc;
}

RVec AtomTree::getInternalForce() {
    calcSpatialForces();
    rbTree.calcTreeInternalForces(spatialForces);
    RVec T(getIVMDim());
    rbTree.getConstraintCorrectedInternalForces(T);
    return T;
}

void AtomTree::enforceConstraints(RVec& pos, RVec& vel) { 
    rbTree.enforceTreeConstraints(pos,vel); 
    rbTree.enforceConstraints(pos,vel);
}

AT_Build::AT_Build( IVM*                ivm,
                    AtomTree*           atree,
                    AtomClusterNode*    acnode,
                    CDSList<AtomLoop>&  aloops)
  : ivm(ivm), loops(aloops), tree(atree), 
    assignedAtoms(tree->ivm->getAtoms().size(),0)
{
    buildAtomClusterNode(acnode);
}

int AtomTree::getIVMDim() const {
    assert(ivm);
    return ivm->dim();
}

// get inclusive (including all children) degrees of freedom + constraints
static int
getIDim(const AtomClusterNode* n) {
    int ret = n->getDim();
    for (int i=0; i < n->getNChildren(); i++)
        ret += getIDim( n->getChild(i) );
    return ret;
}

int 
AtomTree::getDOF() {
    int ret=0;
    for (int i=1 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            ret += nodeTree[i][j]->getDOF();
    return ret;
}

int 
AtomTree::getDim() {
    int              totalDim = 0;
    CDSVector<int,0> moleculeDim;

    if ( nodeTree.size()>1 ) {
        moleculeDim.resize( nodeTree[1].size() );
        for (int i=0 ; i<nodeTree[1].size() ; i++) {
            moleculeDim(i) = getIDim( nodeTree[1][i] );
            totalDim += moleculeDim(i);
        }
    } else
        moleculeDim.resize( 0 );
    return totalDim;
}

//
// Add node to nodeTree. Takes over ownership of the node's heap space.
//
void
AtomTree::addAtomClusterNode(AtomClusterNode* node)
{
    int level = node->level;
    if ( nodeTree.size()<=level )
        nodeTree.resize(level+1);
    nodeTree[level].append(node);
}

class FindBase {   
    IVM*                ivm;
    CDSVector<bool,0>   assignedAtoms; //used to locate cycles while finding base
    IVMAtom*            base;
    FindBase();
    FindBase(const FindBase&);
public:
    ~FindBase() {}
    FindBase& operator=(const FindBase&);

    FindBase(IVM*, IVMAtom* const);
    IVMAtom* recurseToBase(IVMAtom* const,
                           IVMAtom* const);
    operator IVMAtom*() {return base;}
};

FindBase::FindBase(IVM*           ivm,
                   IVMAtom* const i)
  : ivm(ivm), assignedAtoms(ivm->getAtoms().size(),0) 
{
    base = recurseToBase(i,0);
}

//
// Recursive routine to find an end to a tree (a node with only one bond)
// -- for molecules with no atoms fixed in space.
//
// Note that this does not take into account any groupings.
// Side effect: assigned field of atomlist is set to check for cycles.
IVMAtom*
FindBase::recurseToBase(IVMAtom* const a, IVMAtom* const a_prev) {
    assignedAtoms(a->index)=1;
    if (a->bonds.size() == 0) return a; // isolated atom
    if (a->bonds.size() == 1) return a; // a branchless node                
    else {                              // a node with branches
        IVMAtom* a1 = a->bonds[0];
        if (a1==a_prev) 
            a1 = a->bonds[1];
        if ( assignedAtoms(a1->index) ) 
            return a;                   // found a cycle
        return recurseToBase(a1,a);
    }
}

void
AtomTree::addMolecule(AtomClusterNode* groundNode,
                      IVMAtom*         atom ) 
{
    AtomClusterNode* baseNode = new AtomClusterNode(ivm,
                                                    FindBase(ivm,atom),
                                                    groundNode->atoms[0],
                                                    groundNode);
    //FIX: does this value of parentAtom work for torsion hnodes?
    addAtomClusterNode(baseNode);
    groundNode->addChild(baseNode);

    // build up non-fixed proto- hnodes
    AT_Build b(ivm,this,baseNode,loops); 
}

static void
mergeGroups(CDSList< CDSList<int> >& gList) {
    for (int i=0; i<gList.size(); i++)
        for (int k=0; k<gList[i].size(); k++)
            for (int j=i+1; j<gList.size(); j++) {
                if (!gList[j].contains(gList[i][k])) 
                    continue;
                for (int l=0; l<gList[j].size(); l++)
                    if (!gList[i].contains(gList[j][l]))
                        gList[i].append(gList[j][l]);
                gList.remove(j); j--;
            }
}

//
// Construct AtomClusterNode tree consisting of all molecules.
//
AtomTree::AtomTree(IVM* ivm_)
  : ivm(ivm_)
{
    if ( ivm->atoms.size() == 0 ) return;

    CDSVector<bool,0> assignedAtoms( ivm->getAtoms().size() , false );

    loops.resize(0);
    for (int i=0 ; i<ivm->constraintList.size() ; i++)
        loops.append( AtomLoop(ivm->atoms[ ivm->constraintList[i].a ] ,
                               ivm->atoms[ ivm->constraintList[i].b ]) );

    mergeGroups( ivm->groupList );

    //build proto-AtomClusterNodes attached to fixed atoms

    AtomClusterNode* groundNode = new AtomClusterNode( ivm, ivm->atoms[0], 0, 0 );
    AT_Build b(ivm,this,groundNode,loops); 
    addAtomClusterNode(groundNode);
    markAtoms( assignedAtoms );

    //first try trees based on previous base atoms
    for (l_int i=0 ; i<ivm->oldBaseAtoms.size() ; i++)
        if ( !assignedAtoms( ivm->oldBaseAtoms[i] ) ) {
            addMolecule(groundNode,ivm->atoms[ivm->oldBaseAtoms[i]]);
            markAtoms( assignedAtoms );
        }

    //then cover any atoms which are not accounted for
    for (int anum=1 ; anum<ivm->atoms.size() ; anum++)
        // search for an unassigned atom number
        if ( !assignedAtoms( anum ) ) {
            addMolecule(groundNode,ivm->atoms[anum]);
            markAtoms( assignedAtoms );
        }

    ivm->oldBaseAtoms.resize(0);  //remember this set of base atoms
    if ( nodeTree.size()>1 )
        for (int i=0 ; i<nodeTree[1].size() ; i++)
            ivm->oldBaseAtoms.append( nodeTree[1][i]->atoms[0]->index );

    // construct true node data

    int cnt=1; //offset into pos, vel, acc
    nodeTree[0][0] = AtomClusterNode::constructFromPrototype( nodeTree[0][0] , "ground" , cnt);
    // i==1 are base nodes
    for (int i=1 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++) {
            AtomClusterNode* n = nodeTree[i][j];
            HingeSpec hingeSpec("unknown");
            for (int k=0 ; k<ivm->hingeList.size() ; k++)
                //FIX: should check for redundant hinge specifications
                if ( ivm->hingeList[k].aList.contains( n->atoms[0]->index ) )
                    hingeSpec = ivm->hingeList[k];
            nodeTree[i][j] = AtomClusterNode::constructFromPrototype( nodeTree[i][j] , hingeSpec, cnt);
        }

    // Calculate the fixed stations for all the atoms.
    for (int i=0; i<nodeTree.size(); i++)
        for (int j=0; j<nodeTree[i].size(); j++) {
            AtomClusterNode* n = nodeTree[i][j];
            const CDSVec3 OB_G = n->atoms[0]->pos; // origin
            for (int a=0; a < n->atoms.size(); ++a)
                n->atoms[a]->station_B = n->atoms[a]->pos - OB_G;
        }

    ivm->dof_ = getDOF();
    createRigidBodyTree();
}

// sherm: 
// We expect that nodeTree and loops data structures have been filled
// in. Now we create the corresponding rigid body tree and sets of
// length constraints to enforce the loops.
// TODO: using a length (distance) constraint to replace a bond is wrong since
// bonds also have directional constraints. For example, if the bond
// would normally be modeled as a 1-dof torsion joint, their should
// be *5* constraints, not just the single length constraint. Better is
// to leave the bonds alone, split bodies instead, and reconnect them
// with a "weld", which consists of 6 independent constraints.
void AtomTree::createRigidBodyTree() {
    // Go through all nodes from base to tips.
    int nextStateOffset=1; //offset into pos, vel, acc
    for (int i=0; i<nodeTree.size(); i++)
        for (int j=0; j<nodeTree[i].size(); j++) {
            AtomClusterNode& ac = *nodeTree[i][j];
            RigidBodyNode*   rb = RigidBodyNode::create(
                                        ac.getMassPropertiesInBodyFrame(),
                                        ac.getJointFrameInBodyFrame(),
                                        ac.getJointType(), 
                                        ac.getJointIsReversed(),
                                        ivm->minimization(),  // TODO: kludge
                                        nextStateOffset); 
            int rbIndex = -1;
            if (ac.getParent()) {
                RigidBodyNode& parent = rbTree.updRigidBodyNode(ac.getParent()->getRBIndex());
                rbIndex = rbTree.addRigidBodyNode(parent,ac.getReferenceBodyFrameInParent(),rb);
            } else
                rbIndex = rbTree.addGroundNode(rb);

            ac.setRBIndex(rbIndex);
        }

    for (int i=0; i<loops.size(); ++i) {
        AtomLoop& al = loops[i];
        IVMAtom&  t1 = *al.getTip1();
        IVMAtom&  t2 = *al.getTip2();

        RBStation s1(rbTree.updRigidBodyNode(t1.node->getRBIndex()), t1.station_B);
        RBStation s2(rbTree.updRigidBodyNode(t2.node->getRBIndex()), t2.station_B);
        double    d = sqrt(abs2(t2.pos - t1.pos));
        al.setRBDistanceConstraintIndex(rbTree.addDistanceConstraint(s1,s2,d));
    }

    rbTree.finishConstruction(ivm->Ctolerance(), ivm->verbose());
}

void
AtomTree::addCM(const AtomClusterNode* n,
                double&                mass,
                CDSVec3&                  pos)
{
    for (int i=0 ; i<n->atoms.size() ; i++) {
        if ( n->atoms[i]->mass>0.0 ) {
            //   mass += n->atoms[i]->mass;
            mass += 1.0;    //FIX: geometrical center
            //   pos += n->atoms[i]->mass * n->atoms[i]->pos;
            pos += n->atoms[i]->pos;
        }
    }
    for (l_int i=0 ; i<n->children.size() ; i++)
        addCM(n->children[i],mass,pos);
}

//
// Calls addCM to recursively determine center of mass of node and children.
//
CDSVec3
AtomTree::findCM(const AtomClusterNode* n) {
    double mass=0.;
    CDSVec3 pos(0.,0.,0.);
    addCM(n,mass,pos);
    pos /= mass;
    return pos;
}

//
// Mark as assigned all atoms in tree.
//
void
AtomTree::markAtoms(CDSVector<bool,0>& assignedAtoms) {
    assignedAtoms.resize( ivm->atoms.size() );
    for (int j=0 ; j<nodeTree.size() ; j++)
        for (int k=0 ; k<nodeTree[j].size() ; k++)
            for (int l=0 ; l<nodeTree[j][k]->atoms.size() ; l++)
                assignedAtoms( nodeTree[j][k]->atoms[l]->index ) = true;
}


// 
// This routine expects to get a fresh node (body) with only one 
// atom attached (the atom at which the inboard joint attaches).
// We'll look through all the defined groups to see if there is one
// containing our atom, and if so we will assign all the other atoms
// in that group to this node.
//
// In the current implementation 
// each atom is allowed to be a member of a single group. All
// members of the group are held rigid wrt each other.
//
void
AT_Build::buildAtomClusterNode(AtomClusterNode* node)
{
    assert(node->atoms.size()==1);

    for (int i=0 ; i<ivm->groupList.size() ; i++)
        if ( ivm->groupList[i].contains( node->atoms[0]->index ) )
            for (int j=0 ; j<ivm->groupList[i].size() ; j++) 
                if ( 0 <= ivm->groupList[i][j] && ivm->groupList[i][j] <  ivm->atoms.size() ) {
                    IVMAtom* atom = ivm->atoms[ ivm->groupList[i][j] ];
                    if ( !node->atoms.contains(atom) ) {
                        node->atoms.append( atom );
                        atom->node = node;
                    }
                }

    // note that all atoms in this node have been assigned
    for (l_int i=0 ; i<node->atoms.size() ; i++)
        assignedAtoms( node->atoms[i]->index ) = true;

    int numAtoms = node->atoms.size(); // we might add atoms- don't process twice
    for (l_int i=0 ; i<numAtoms ; i++) {
        IVMAtom* atom = node->atoms[i];

        for (int j=0 ; j<atom->bonds.size() ; j++) {
            IVMAtom* cAtom = atom->bonds[j];
            // Skip bonded atoms within this node.
            if (cAtom != node->parentAtom && !node->atoms.contains(cAtom)) {
                if ( assignedAtoms(cAtom->index) ) {
                    // cycle - deal with it
                    //if ( InternalDynamics::verbose&InternalDynamics::printLoopInfo )
                    cout << "AT_Build::buildNode: cycle link found between atoms " 
                         << ivm->idAtom(atom->index) << " and " 
                         << ivm->idAtom(cAtom->index) << '\n'
                         << "\tremoving bond." << endl;

                    if ( ivm->useLengthConstraints() ) {
                        loops.append( AtomLoop(cAtom,atom) );
                        cout << "\tadding constraint\n" ;
                    }

                    if ( cAtom->bonds.contains(atom) )
                        cAtom->bonds.remove( cAtom->bonds.getIndex(atom) ) ;
                    else
                        cout << "AT_Build::buildNode: cycle link: unable to remove 2nd bond"
                             << "\n";

                    atom->bonds.remove(j);j--;
                } else {
                    AtomClusterNode* pNode = node;
                    IVMAtom*         pAtom = atom;
                    // bend and bendtorsion hinges require splitting an atom into two.
                    for (int k=0 ; k < ivm->hingeList.size() ; k++)
                        if ((ivm->hingeList[k].type == "bendtorsion" || ivm->hingeList[k].type == "bend")
                            && ivm->hingeList[k].aList.contains(cAtom->index) && pAtom->pos != cAtom->pos )
                        {
                            // check that bend atoms are bound
                            if (   !cAtom->bonds.contains(ivm->atoms[ivm->hingeList[k].atom0])
                                || !cAtom->bonds.contains(ivm->atoms[ivm->hingeList[k].atom1]) )
                            {
                                CDSString mesg = CDSString("AT_Build::buildNode: ") +
                                                    ivm->hingeList[k].type +
                                                    ": specified atoms not bound.";
                                cerr << mesg << '\n'
                                    << "\thinge atom: " << atom << endl
                                    << "\tother atoms: " 
                                    << ivm->atoms[ivm->hingeList[k].atom0]
                                    << " " << ivm->atoms[ivm->hingeList[k].atom1] << endl;
                                throw Exception(mesg);
                            }

                            // pAtom -- cAtom  --> pAtom -- (newAtom -- cAtom)
                            IVMAtom* newAtom = new IVMAtom(ivm->atoms.size(), 0.5*cAtom->mass);
                            cAtom->mass *= 0.5;
                            ivm->atoms.append( newAtom );
                            assignedAtoms.resize( ivm->atoms.size() );  //we've counted it.
                            assignedAtoms( ivm->atoms.size()-1 ) = true;

                            // make a couple of new bonds.
                            newAtom->bonds.append(pAtom);     
                            pAtom->bonds[j] = newAtom; //replace this bond
                            newAtom->bonds.append(cAtom);
                            cAtom->bonds.append(newAtom);

                            //delete the old bond
                            cAtom->bonds.remove( cAtom->bonds.getIndex(pAtom) );
                            //pAtom->bonds.remove( pAtom->bonds.getIndex(cAtom) ); ???

                            newAtom->pos = cAtom->pos;
                            newAtom->vel = cAtom->vel;
                            // bendtorsion hinges require insertion of a new node containing 
                            // one of the resulting atoms
                            if (ivm->hingeList[k].type == "bendtorsion") {
                                InternalDynamics::HingeSpec hingeSpec;
                                hingeSpec.type = "torsion";
                                hingeSpec.aList.append( newAtom->index );
                                ivm->hingeList.append( hingeSpec );

                                AtomClusterNode* newNode = new AtomClusterNode(ivm, newAtom, pAtom, pNode);
                                tree->addAtomClusterNode(newNode);
                                pNode->addChild(newNode);
                                pNode = newNode;
                            } else if (ivm->hingeList[k].type == "bend") {
                                pNode->atoms.append(newAtom);
                            }

                            pAtom = newAtom;
                            break;
                        }

                    // the local node and current atom are remote to newNode
                    AtomClusterNode* newNode = new AtomClusterNode(ivm, cAtom, pAtom, pNode);
                    tree->addAtomClusterNode(newNode);
                    pNode->addChild(newNode);
                    //     assignedAtoms(cAtom->index) = 1;
                    buildAtomClusterNode(newNode);
                }
            }
        }
    }
}

ostream& 
operator<<(ostream& s, const AtomTree& aTree) {
    for (int i=0 ; i<aTree.nodeTree.size() ; i++)
        for (int j=0 ; j<aTree.nodeTree[i].size() ; j++) {
            const AtomClusterNode& ac = *aTree.nodeTree[i][j];
            const RigidBodyNode& rb = aTree.rbTree.getRigidBodyNode(ac.getRBIndex());
            s << "\tnode " << i << ' ' << j 
            << ": " << ac << ' ' << ac.type() << '\n';
            //rb.nodeDump(s);
        }
    return s;
}

typedef CDSVector<CDSVec3> VecVec3;


//  (see CDS notes of 3/27/00)
//
// 1) calculate ``spatial velocities'' from atom velocities:
//   w_i = sum_j cross(q_ij - q_i Va_ij) 
//   v_i = sum_j  Va_ij 
// 2) calculate ``internal velocities'' from these velocites using
//    Vi_i = H_i( Vs_i -  phi_(i-1,i) Vs_i-1 )
// 3) set the internal force T to his value.
// 4) Build P, replacing the mass M_i matrix with R_iR_i^T
//    where R_i = [ 0 1 ... -(q_ij-q_i) 1]
// 5) set the cartesian force, a and b to zero.
// 6) the resulting calculated acceleration is used at the internal 
//    velocities.
//
//
void
AtomTree::velFromCartesian(const RVec& pos, RVec& vel)
{
    VecVec3 avel0(ivm->atoms.size()-1);   // save desired atom velocities
    for (int i=1 ; i<ivm->atoms.size() ; i++)
        avel0(i-1) = ivm->atoms[i]->vel;

    setPos(pos);        // set configuration and calculate related kinematics (incl. atom pos)
    vel.set(0.0);
    rbTree.setVel(vel); // zero velocities to nuke bias forces (no effect on atom vel)

    // calculate impulses from desired atomic momenta and convert to internal coordinates
    calcP();
    calcSpatialImpulses(); // sets spatialForces vector
    rbTree.calcZ(spatialForces);

    // solve for desired internal velocities
    vel = calcGetAccel();
    setVel(vel);    // sets atom vels too
    fixVel0(vel);

    if ( ivm->verbose()&printVelFromCartCost ) {
        VecVec3 avel(ivm->atoms.size()-1);
        for (l_int i=1 ; i<ivm->atoms.size() ; i++)
            avel(i-1) = ivm->atoms[i]->vel;

        // calculate cost
        const double cost = abs2(avel - avel0) / avel.size();
        cout << "velFromCartesian: cost: " << cost << '\n';
    }
}
