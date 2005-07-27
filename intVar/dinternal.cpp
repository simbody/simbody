/**@file
 *
 *  Dynamics in internal coordinates
 * 
 *  ****add references****
 * 
 *  8/25/99 - C. Schwieters
 * 
 *  outline:
 *  
 *  -need simple constraint mechanism (for these coordinates)
 *  -add more general constraint mechanism (for arbitrary constraints)
 *  -convert from internal to cartesian coordinates
 */

#include "dinternal.h"

#include "dint-node.h"
#include "dint-step.h"
#include "dint-loop.h"

#include "linemin.h"
#include "dint-atom.h"
#include "vec3.h"
#include "vec4.h"

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


#include <cdsAuto_arr.h>

using namespace InternalDynamics;
using MatrixTools::inverse;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat6;


//
// structure of atomtree:
//  doubly linked tree
//

class AT_Build {  //temporary class used in building molecule trees
public:
    IVM* ivm;
    CDSList<Loop>& loops;
    AtomTree* const tree;
    CDSVector<bool,0> assignedAtoms; //used in locating cycle links
    AT_Build(IVM*,
             AtomTree* const, 
             HingeNode*,
             CDSList<Loop>&);
    void buildNode(HingeNode*);
};

AtomTree::~AtomTree() {
    for (int i=0 ; i<nodeTree.size() ; i++) {
        for (int j=0 ; j<nodeTree[i].size() ; j++) 
            delete nodeTree[i][j];
        nodeTree[i].resize(0);
    }
    nodeTree.resize(0);
}


//
// destructNode - should also work on partially constructed objects
//
void
AtomTree::destructNode(HingeNode *n) 
{
    for (int i=0 ; i<n->children.size() ; i++)
        destructNode( n->children[i] );
    nodeTree[n->level].remove( nodeTree[n->level].getIndex(n) );
    delete n;
}

AT_Build::AT_Build( IVM*            ivm,
                    AtomTree* const tree,
                    HingeNode*      node,
                    CDSList<Loop>&  loops)
  : ivm(ivm), loops(loops), tree(tree), 
    assignedAtoms(tree->ivm->getAtoms().size(),0)
{
    buildNode(node);
}

// get inclusive (including all children) degrees of freedom + constraints
static int
getIDim(const HingeNode* n) {
    int ret = n->getDim();
    for (int i=0 ; n->getChild(i) ; i++)
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
    int ret = 0;
    CDSVector<int,0> moleculeDim;

    if ( nodeTree.size()>1 ) {
        moleculeDim.resize( nodeTree[1].size() );
        for (int i=0 ; i<nodeTree[1].size() ; i++) {
            moleculeDim(i) = getIDim( nodeTree[1][i] );
            ret += moleculeDim(i);
        }
    } else
        moleculeDim.resize( 0 );
    return ret;
}

//
// Add node to nodeTree.
//
void
AtomTree::addNode(HingeNode* node)
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
AtomTree::addMolecule(HingeNode* originNode,
                      IVMAtom*   atom ) 
{
    HingeNode* baseNode = new HingeNode(ivm,
                                        FindBase(ivm,atom),
                                        originNode->atoms[0],
                                        originNode);
    //FIX: does this value of remAtom work for torsion hnodes?
    addNode(baseNode);
    originNode->addChild(baseNode);

    // build up non-fixed proto- hnodes
    AT_Build b(ivm,this,baseNode,ivm->loops); 
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
// Construct tree consisting of all molecules.
//
AtomTree::AtomTree(IVM* ivm) : ivm(ivm)
{
    if ( ivm->atoms.size() == 0 ) return;

    CDSVector<bool,0> assignedAtoms( ivm->getAtoms().size() , 0 );
    ivm->loops.resize(0);
    for (int i=0 ; i<ivm->constraintList.size() ; i++)
        ivm->loops.append( Loop(ivm->atoms[ ivm->constraintList[i].a ] ,
                                ivm->atoms[ ivm->constraintList[i].b ]) );

    mergeGroups( ivm->groupList );

    //build proto-hnodes attached to fixed atoms

    HingeNode* originNode = new HingeNode( ivm, ivm->atoms[0] );
    AT_Build b(ivm,this,originNode,ivm->loops); 
    addNode(originNode);
    markAtoms( assignedAtoms );

    //first try trees based on previous base atoms
    for (l_int i=0 ; i<ivm->oldBaseAtoms.size() ; i++)
        if ( !assignedAtoms( ivm->oldBaseAtoms[i] ) ) {
            addMolecule(originNode,ivm->atoms[ivm->oldBaseAtoms[i]]);
            markAtoms( assignedAtoms );
        }

    //then cover any atoms which are not accounted for
    for (int anum=1 ; anum<ivm->atoms.size() ; anum++)
        // search for an unassigned atom number
        if ( !assignedAtoms( anum ) ) {
            addMolecule(originNode,ivm->atoms[anum]);
            markAtoms( assignedAtoms );
        }

    ivm->oldBaseAtoms.resize(0);  //remember this set of base atoms
    if ( nodeTree.size()>1 )
        for (int i=0 ; i<nodeTree[1].size() ; i++)
            ivm->oldBaseAtoms.append( nodeTree[1][i]->atoms[0]->index );

    // construct true node data

    int cnt=1; //offset into pos, vel, acc
    nodeTree[0][0] = construct( nodeTree[0][0] , "origin" , cnt);
    // i==1 are base nodes
    for (l_int i=1 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++) {
            HingeNode* n = nodeTree[i][j];
            HingeSpec hingeSpec("unknown");
            for (int k=0 ; k<ivm->hingeList.size() ; k++)
                //FIX: should check for redundant hinge specifications
                if ( ivm->hingeList[k].aList.contains( n->atoms[0]->index ) )
                    hingeSpec = ivm->hingeList[k];
            nodeTree[i][j] = construct( nodeTree[i][j] , hingeSpec, cnt);
        }

    ivm->dof_ = getDOF();
    ivm->lConstraints->construct(ivm->loops);
}

void
AtomTree::addCM(const HingeNode* n,
                double&          mass,
                Vec3&            pos)
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
Vec3
AtomTree::findCM(const HingeNode* n) {
    double mass=0.;
    Vec3 pos(0.,0.,0.);
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
                assignedAtoms( nodeTree[j][k]->atoms[l]->index ) = 1;
}


//
// Add other atoms from group here. In the current implementation 
// each atom is allowed to be a member of a single group. All
// members of the group are held rigid wrt each other.
//
void
AT_Build::buildNode(HingeNode* node)
{
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

    for (l_int i=0 ; i<node->atoms.size() ; i++)
        assignedAtoms( node->atoms[i]->index ) = 1;

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
                        loops.append( Loop(cAtom,atom) );
                        cout << "\tadding constraint\n" ;
                    }
                    if ( cAtom->bonds.contains(atom) )
                        cAtom->bonds.remove( cAtom->bonds.getIndex(atom) ) ;
                    else
                        cout << "AT_Build::buildNode: cycle link: unable to remove 2nd bond"
                            << "\n";
                    atom->bonds.remove(j);j--;
                } else {
                    HingeNode* pNode = node;
                    IVMAtom*   pAtom = atom;
                    // bend and bendtorsion hinges require splitting an atom into two.
                    for (int k=0 ; k < ivm->hingeList.size() ; k++)
                        if ((ivm->hingeList[k].type == "bendtorsion" || ivm->hingeList[k].type == "bend")
                            && ivm->hingeList[k].aList.contains(cAtom->index) && pAtom->pos != cAtom->pos )
                        {
                            // check that bend atoms are bound
                            if (   !cAtom->bonds.contains(ivm->atoms[ivm->hingeList[k].atom0])
                                || !cAtom->bonds.contains(ivm->atoms[ivm->hingeList[k].atom1]) )
                            {
                                String mesg = String("AT_Build::buildNode: ") +
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
                            assignedAtoms( ivm->atoms.size()-1 ) = 1;

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

                                HingeNode* newNode = new HingeNode(ivm, newAtom, pAtom, pNode);
                                tree->addNode(newNode);
                                pNode->addChild(newNode);
                                pNode = newNode;
                            } else if (ivm->hingeList[k].type == "bend") {
                                pNode->atoms.append(newAtom);
                            }

                            pAtom = newAtom;
                            break;
                        }

                    // the local node and current atom are remote to newNode
                    HingeNode* newNode = new HingeNode(ivm, cAtom, pAtom, pNode);
                    tree->addNode(newNode);
                    pNode->addChild(newNode);
                    //     assignedAtoms(cAtom->index) = 1;
                    buildNode(newNode);
                }
            }
        }
    }
}

ostream& 
operator<<(ostream& s, const AtomTree& aTree) {
    for (int i=0 ; i<aTree.nodeTree.size() ; i++)
        for (int j=0 ; j<aTree.nodeTree[i].size() ; j++)
            s << "\tnode " << i << ' ' << j 
            << ": " << *aTree.nodeTree[i][j] << ' '
            << aTree.nodeTree[i][j]->type() << '\n';
    return s;
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
AtomTree::calcP() {
    // level 0 for atoms whose position is fixed
    // InternalDynamics::RecurseTipToBase<HingeNode::calcPandZ> d1;
    for (int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcP();
}


// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
AtomTree::calcZ() {
    // level 0 for atoms whose position is fixed
    for (int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcZ();
}


// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
AtomTree::calcPandZ() {
    // level 0 for atoms whose position is fixed
    // InternalDynamics::RecurseTipToBase<HingeNode::calcPandZ> d1;
    for (int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcPandZ();
}


IVM::IVM() 
  : tree_(0),                  integrator_(0),
    dof_(0),                   dim_(0),            verbose_( printStepInfo ), 
    useLengthConstraints_(0),  adjustTS_(1),       scaleVel_(1),
    resetCMflag(0),
    Etotal_(0.0),              Epotential_(0.0),   Ekinetic_(0.0),
    Etolerance_(0),            Gtolerance_(1e-6),  Ctolerance_(1e-8),
    minStepSize_(1e-10),       maxTSFactor_(1.025),
    maxDeltaE_(100),           maxCalls_(5000),
    kBoltzmann_(0.0),          bathTemp_(298.0),   frictionCoeff_(0.0),
    responseTime_(0.0),        currentTemp_(-1.0),
    dEpred_(0.001),            integrateType( "PC6" ),
    rvecSize_(defaultRVecSize), rvecProd_(defaultRVecProd),
    vecVec3Size_(defaultVecVec3Size), vecVec3Prod_(defaultVecVec3Prod)
{
    lConstraints = new LengthConstraints(this);
    tree_ = new AtomTree(this);
}

IVM::~IVM()
{
    delete lConstraints;
    delete integrator_;
    delete tree_;
    for (int i=0 ; i<atoms.size() ; i++)
        delete atoms[i];
}

//
// Y is used for length constraints.
//
void
IVM::calcY() {
    for (int i=0 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++)
            tree()->nodeTree[i][j]->calcY();
}

bool
IVM::minimization() const { 
    if ( !integrator_ ) {
        cout << "IVM::minimization: integrator not yet defined.\n" << flush;
        throw Exception("IVM::minimization: integrator not yet defined.");
    }
    return integrator()->minimization(); 
}

// Calc acceleration: sweep from base to tip.
RVec
AtomTree::calcGetAccel() {
    RVec acc( ivm->dim() );
    for (int i=0 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcAccel();
    for (l_int i=0 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->getAccel(acc);
    return acc;
}

// Calc P quantities: sweep from tip to base.
RVec
AtomTree::getAccel() {
    calcPandZ();
    RVec acc = calcGetAccel();

    if ( ivm->lConstraints->fixAccel() ) 
        for (int i=0 ; i<nodeTree.size() ; i++)
            for (int j=0 ; j<nodeTree[i].size() ; j++)
                nodeTree[i][j]->getAccel(acc);

    return acc;
}

void
AtomTree::updateAccel() {
    for (int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcZ();
    for (l_int i=0 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcAccel();
}

// Calc internal force: sweep from tip to base.
RVec
AtomTree::getInternalForce() {
    for (int i=nodeTree.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->calcInternalForce();

    RVec T( ivm->dim() );
    for (l_int i=0 ; i<nodeTree.size() ; i++)
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->getInternalForce(T);

    ivm->lConstraints->fixGradient(T);

    return T;
}

void 
AtomTree::setPosVel(const RVec& pos, const RVec& vel) {
    for (int i=0 ; i<nodeTree.size() ; i++) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->setPosVel(pos,vel); 
}

void 
AtomTree::setVel(const RVec& vel)  {
    for (int i=0 ; i<nodeTree.size() ; i++) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->setVel(vel); 
}

void 
AtomTree::enforceConstraints(RVec& pos, RVec& vel) {
    for (int i=0 ; i<nodeTree.size() ; i++) 
        for (int j=0 ; j<nodeTree[i].size() ; j++) 
            nodeTree[i][j]->enforceConstraints(pos,vel);

    ivm->lConstraints->enforce(pos,vel); //FIX: previous constraints still obeyed?
}

RVec
AtomTree::getPos() const {
    RVec pos(ivm->dim());
    for (int i=0 ; i<nodeTree.size() ; i++) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->getPos(pos); 
    return pos;
}

RVec
AtomTree::getVel() const {
    RVec vel(ivm->dim());
    for (int i=0 ; i<nodeTree.size() ; i++) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->getVel(vel); 
    return vel;
}

int
IVM::defaultRVecSize(const RVec& v) {
    return v.size();
}

double
IVM::defaultRVecProd(const RVec& v1, const RVec& v2) {
    return dot(v1,v2);
}

int
IVM::defaultVecVec3Size(const CDSVector<Vec3>& v) {
    return v.size(); // FIX: *3???
}

double
IVM::defaultVecVec3Prod(const CDSVector<Vec3>& v1, const CDSVector<Vec3>& v2) {
    // return dot(v1,v2);// didn't work??
    double sum=0.;
    for (int i=0 ; i<v1.size() ; i++)
        sum += v1[i]*v2[i];
    return sum;
}

//
// recursively group all atoms in a node with neighboring nodes 
// appropriate for torsion angle dynamics. 
// does not create a group if there's only a single atom.
//
void
IVM::groupTorsion(const HingeNode* n) {
    CDSList<int> aList;
    const IVMAtom* a;

    //fixed atoms
    if ( n->isOriginNode() ) {
        if ( n->getAtom(1) ) 
            for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                if ( a->index>=0 )
                    aList.append( a->index );

        const HingeNode* c;
        for (int i=0 ; ( c=n->getChild(i) ) ; i++) { 
            if ( c->getChild(0) || c->getAtom(1) )     
                //has children or contains more than one atom
                groupTorsion( c );       
            else
                if (c->parentAtom->index!=0 && c->getAtom(0)->index>=0)
                    //one atom, no children and not free atom : group
                    aList.append( c->getAtom(0)->index );
        }
    } else { //not origin node
        //base node w/ 1 or 2 atoms- w/ no bond to base and one or fewer children
        if ( n->isBaseNode() && !n->getAtom(2) &&
             n->parentAtom->index==0 && !n->getChild(1) )
        { 
            if ( n->getAtom(1) || n->getChild(0) )
                for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                    if ( a->index>=0 )
                        aList.append( a->index );
            n = n->getChild(0);
        }

        if (n) {
            const HingeNode* c;
            for (int i=0 ; ( c=n->getChild(i) ) ; i++) {
                if ( c->getChild(0) || c->getAtom(1) ) {    
                    //has children or contains more than one atom
                    groupTorsion( c );       
                } else {
                    //no children: group
                    for (int j=0 ; ( a=c->getAtom(j) ) ; j++)
                        if ( a->index>=0 ) {
                            bool add = 1;
                            for (int k=0 ; k<hingeList.size() ; k++)
                                if ( hingeList[k].aList.contains(a->index) &&
                                     hingeList[k].type != "torsion")
                                {
                                    add = 0;
                                    break;
                                }
                            if ( add )
                                aList.append( a->index );
                        }
                }
            }

            if ( aList.size() || n->getAtom(1) )
                for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                    if ( a->index>=0 )
                        aList.append( a->index );
        }
    }

    if ( aList.size() ) groupList.append( aList );
}

//
// GroupList is overwritten, but pre-existing groupings are honored.
//
void
IVM::groupTorsion()
{
    groupList.resize(0);   //clobber existing groupList
    groupTorsion( tree()->nodeTree[0][0] );
}

void
IVM::step(double& timeStep) {
    bool ok=0;
    do {
        try {
            integrator()->step(timeStep);
            ok=1;
        }
        catch ( Exception a ) {
            if ( String(a.mess).contains("large time-step") ) {
                cout << "InternalDynamics::step: large timestep detected. Halving.\n";
                timeStep *= 0.5;
                if ( timeStep < minStepSize() ) {
                    cout << "InternalDynamics::step: stepsize too small." << endl;
                    throw InternalDynamics::Exception("stepsize too small");
                }
            } else {
                cout << "InternalDynamics::step: caught assertion failure:\n"
                << "\t" << a.mess << '\n';
                cout.flush();
                throw;
            }
        }
    } while ( !ok );
}
 
//void
//InternalDynamics::approxTemperature()
//{
// approxTemp = 0.0;
// for (int i=0 ; i<tree->nodeTree.size() ; i++)
//   for (int j=0 ; j<tree->nodeTree[i].size() ; j++) 
//     approxTemp += tree->nodeTree[i][j]->approxKE();
//     
// approxTemp /= (dof * kBoltzman);
//} /* calcTemp */
 
void
IVM::calcTemperature() {
    Ekinetic_ = 0.0;
    int enumerateKE=0;
    if (enumerateKE) 
        cout << "thermal energy contributions (by node):\n";

    for (int i=0 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
            if (enumerateKE) 
                cout << i << ' ' << j << ' ';
            double ke = tree()->nodeTree[i][j]->kineticE();
            if (enumerateKE) 
                cout << ke << '\n';
            Ekinetic_ += ke;
        }


    if ( dof() )
        currentTemp_ = 2.0 * Ekinetic() / (dof() * kBoltzmann());
    else
        currentTemp_ = 0;
    // approxTemperature(); //FIX: temporary. remove after testing.
}
 
typedef SubVector<Vec6> RSubVec6;
typedef SubVector<const Vec6> ConstRSubVec6;

//
// Using current atom positions/velocities.
// This resets the CM vel and angular momentum. This operation is performed 
// in Cartesian space, and then the internal coordinates are recalculated.
// Thus, the linear vel. and AM will not generally be exactly zero after
// converting to internal coords and reconverting to Cartesians.
//
// do this in spatial coords. L = L_l + L_a
//
// correct base velocity and then recalc spatial and cartesian vels.
//
// Remove only the linear velocity component from each body. 
void
IVM::resetCM()
{
    if ( !resetCMflag ) return;
    resetCMflag=0;

    // exit if there are any fixed atoms: in this case, 
    //  removing CM velocity will change relative (internal) velocities
    if ( tree()->nodeTree[0][0]->getAtom(1) ) return;

    static AtomList nodes; //to contain mass, posCM, velCM of each node
    nodes.resize(0);
    for (int i=1 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
            const HingeNode* n = tree()->nodeTree[i][j];
            IVMAtom *a = new IVMAtom( nodes.size() , n->mass() );
            a->pos = n->posCM();
            a->vel = Vec3( ConstRSubVec6(n->getSpatialVel(),3,3).vector() )
                     - cross( ConstRSubVec6(n->getSpatialVel(),0,3).vector() , 
                              n->getAtom(0)->pos - n->posCM() ); 
            nodes.append( a );
        }

    double mass=0.0;             //overall mass, posCM, velCM
    Vec3 posCM(0.0), velCM(0.0);
    for (l_int i=0 ; i<nodes.size() ; i++) {
        mass += nodes[i]->mass;
        posCM += nodes[i]->mass * nodes[i]->pos;
        velCM += nodes[i]->mass * nodes[i]->vel;
    }

    posCM /= mass;
    velCM /= mass;

    Vec3 L(0.0); // angular momentum excluding internal rotations
    for (l_int i=0 ; i<nodes.size() ; i++)
        L += nodes[i]->mass * cross( nodes[i]->pos - posCM , 
                                     nodes[i]->vel - velCM );

    if (verbose_&printResetCM) 
        cout << "resetCM: velCM: " << velCM << "overall AM: " << L << '\n';
    Vec3 omega(0.0); // the overall rotational velocity
    if ( nodes.size() >2 ) {
        InertiaTensor inertia;
        inertia.calc( posCM , nodes );
        try {
            omega = inverse(inertia) * L;
        } 
        catch ( CDS_NAMESPACE(SingularError) ) {
            cout << "InternalDynamics::resetCM: singular inertia Tensor.\n";
        }
    } else if ( nodes.size() == 2) {
        Vec3 r = nodes[1]->pos - nodes[0]->pos;
        omega = 1.0/abs2(r) * cross( r , nodes[1]->vel - nodes[0]->vel );
    }   

    // FIX: if a node lies on a principle axis
    //   add the appropriate component(s) to the total AM;

    // integrator()->setPos( tree()->getPos() ); //FIX: ??
    // integrator()->setVel( tree()->getVel() );

    for ( l_int i=0 ; i<tree()->nodeTree[1].size() ; i++ ) {
        HingeNode* n = tree()->nodeTree[1][i];
        Vec6 sVel = n->getSpatialVel() - blockVec(omega, 
                                          cross( omega , 
                                                 n->getAtom(0)->pos-posCM) + velCM);
        n->setVelFromSVel(sVel);
    }
    integrator()->setPos( tree()->getPos() );
    integrator()->setVel( tree()->getVel() );

    tree()->enforceConstraints(integrator()->getPos(),integrator()->getVel());
    tree()->setPosVel(integrator()->getPos(),integrator()->getVel());

    if (verbose_&printResetCM) { //REMOVE: for testing only
        Vec3 velCM(0.0);
        for (l_int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
                const HingeNode* n = tree()->nodeTree[i][j];
                Vec3 vCM = Vec3( ConstRSubVec6(n->getSpatialVel(),3,3).vector() )
                                 - cross( ConstRSubVec6(n->getSpatialVel(),0,3).vector() , 
                                          n->getAtom(0)->pos - n->posCM() ); 
                velCM += n->mass() * vCM;
            }
        velCM /= mass;

        Vec3 L(0.0); // angular momentum excluding internal rotations
        for (l_int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
                const HingeNode* n = tree()->nodeTree[i][j];
                Vec3 vCM = Vec3( ConstRSubVec6(n->getSpatialVel(),3,3).vector() ) - 
                                cross( ConstRSubVec6(n->getSpatialVel(),0,3).vector() , 
                                        n->getAtom(0)->pos - n->posCM() ); 
                L += n->mass() * cross( n->posCM() - posCM , vCM - velCM );
            }
        cout << "resetCM: after reset: velCM: " << velCM << "overall AM: " << L << '\n';
    }

    for (l_int i=0 ; i<nodes.size() ; i++) delete nodes[i]; //clean up.

    // --these next two lines are necessary if the pot. depends on absolute
    //   position or orientation
    // calcEnergy();        //recalc these- so eq. of motion are correct
    // calcAccel(*trees);
}

void
IVM::printCM()
{
    if ( verbose_&printCMVel ) {
        Vec3 posCM(0.0);
        double mass=0.0;;
        for (int i=1 ; i<atoms.size() ; i++) {
            IVMAtom* a = atoms[i];
            posCM += a->mass * a->pos;
            mass += a->mass;
        }
        posCM /= mass;

        Vec3 velCM(0.0);
        for (l_int i=1 ; i<atoms.size() ; i++)        //remove CM velocity
            velCM += atoms[i]->mass * atoms[i]->vel;
        velCM /= mass;

        cout.setf( ios::fixed | ios::right);

        cout << "COM vel: " << setw(6) << setprecision(3) << velCM << ' ';

        Vec3 L(0.0);
        for (l_int i=1 ; i<atoms.size() ; i++)
            L += atoms[i]->mass * cross( atoms[i]->pos-posCM , atoms[i]->vel );

        cout << "ang mom: " << setw(6) << L << '\n';
    }

    //FIX: this should be removed
    calcTemperature();

    if ( verbose_&printTemperature )
        cout << "temperature: " << currentTemp() << '\n';
    if ( verbose_&printEnergy )
        cout << "Energy: " << setprecision(14) << Etotal() 
             << " (" << Epotential() << ' ' << Ekinetic() << ')'
             << '\n';
    if ( verbose_&printCoords )
        cout << "pos: " << setprecision(14) << integrator()->getPos() << '\n';
    if ( verbose_&(printNodeForce | printNodePos | printNodeTheta) )
        for (int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++)
                tree()->nodeTree[i][j]->print(verbose_);
}

// FIX: this is really xplor-specific
void
IVM::initAtoms(const int     natom,
               const double* massA,
               const double& kB)
{
    MALLOC_DEBUG;
    INSTALL_FPE_HANDLER;

    kBoltzmann_ = kB;
    for (int i=0 ; i<atoms.size() ; i++)
        delete atoms[i];
    atoms.resize(natom+1);
    atoms[0] = new IVMAtom(0,0.0);   //origin atom
    for (l_int i=1 ; i<=natom ; i++)
        atoms[i] = new IVMAtom(i,massA[i-1]);
    loops.resize(0);
    cout.precision(8);
}

//void
//InternalDynamics::initTree()
//{
//} /* initTrees */

void
IVM::updateAccel() {
    tree()->updateAccel();
}

void
IVM::initTopology()
{
    // integrator must be defined to determine hinge setup
    //  (whether to use euler angles or quaternions)
    if (integrator()) delete integrator_;
    integrator_ = Integrator::create(integrateType,this);

    if (tree())
        delete tree_;
    tree_ = new AtomTree(this);

    if ( verbose_&printNodeDef ) {
        cout << "Node info:\n";
        cout << *tree() << '\n';
        cout << "total number of degrees of freedom: " << dof() << '\n';
    }

    dim_ = tree()->getDim();
}

void
IVM::initDynamics(bool reuseTopology)
{
    RVec pos, vel, acc;
    if ( !reuseTopology ) initTopology();

    pos = tree()->getPos();
    if ( minimization() ) {
        vel.resize( pos.size() );
        vel.set(0.0);
        tree()->enforceConstraints(pos,vel);    //FIX: probably redundant
        tree()->setPosVel(pos,vel);
    } else { // integration
        vel = tree()->getVel();          
        tree()->velFromCartesian(pos,vel);
        tree()->enforceConstraints(pos,vel);    //FIX: probably redundant
        resetCM();
        tree()->setPosVel(pos,vel);
    }

    calcEnergy();
    if ( !minimization() )
        acc = tree()->getAccel();
    integrator()->init(pos,vel,acc);
}

String
IVM::idAtom(int id) const {
    StringStream ret;
    ret << id << ends;
    return ret.str();
}

typedef CDSVector<Vec3> VecVec3;

//
// set up sVel = R * M * Va (CDS notes of 3/27/00)
//
void
AtomTree::propagateSVel()
{
    for (int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++) {
            HingeNode* n = nodeTree[i][j];
            Vec6 sVel(0.0);
            for (int k=0 ; n->getAtom(k) ; k++)
                sVel += n->getAtom(k)->mass
                        * blockVec( cross(n->getAtom(k)->pos - n->getAtom(0)->pos,
                                          n->getAtom(k)->vel) ,
                                    n->getAtom(k)->vel );
            nodeTree[i][j]->propagateSVel( sVel );
        }
}

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
    VecVec3 avel0(ivm->atoms.size()-1);
    for (int i=1 ; i<ivm->atoms.size() ; i++)
        avel0(i-1) = ivm->atoms[i]->vel;

    vel.set(0.0);    //   can we set just pos instead?
    setPosVel(pos, vel);

    for (l_int i=1 ; i<ivm->atoms.size() ; i++) //reset atom velocities
        ivm->atoms[i]->vel = avel0(i-1);

    for (l_int i=nodeTree.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<nodeTree[i].size() ; j++)
            nodeTree[i][j]->prepareVelInternal();

    // turn off Langevin stuff for now
    double frictionCoeff = ivm->frictionCoeff();
    ivm->frictionCoeff_ = 0.0;

    // make sure that internal velocities of each node are correct

    propagateSVel();
    calcPandZ();

    vel = calcGetAccel();
    setVel(vel);
    ivm->lConstraints->fixVel0(vel);

    VecVec3 avel(ivm->atoms.size()-1);
    for (l_int i=1 ; i<ivm->atoms.size() ; i++)
        avel(i-1) = ivm->atoms[i]->vel;

    // calculate cost
    double cost = abs2(avel - avel0) / avel.size();
    if ( ivm->verbose()&printVelFromCartCost )
        cout << "velFromCartesian: cost: " << cost << '\n';

    ivm->frictionCoeff_ = frictionCoeff;
}

//ostream& 
//InternalDynamics::operator<<(      ostream&  s,
//                   const InternalDynamics::HingeSpec& h)
//{
// // long oflags = s.flags();
// // s.setf( ios::left );
// s << " " << setw(9) << h.type << "   " << h.aList << ' ';
// if ( h.atom0>=0 )
//   s << "( " << h.atom0 << " ) ";
// if ( h.atom1>=0 )
//   s << "( " << h.atom1 << " )";
// // s << '\n';
// // s.flags( oflags );
// return s;
//} /* operator<<(HingeSpec) */
//
