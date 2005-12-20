/**@file
 * This is the implementation of the AtomClusterNode class used in constructing
 * a rigid body model from a set of atoms, bonds, and miscellaneous instructions.
 */

#include "dinternal.h"
#include "AtomTree.h"
#include "AtomClusterNode.h"

#include "dint-step.h"
#include "dint-atom.h"

#include "cdsMath.h"
#include "cdsVector.h"
#include "fixedVector.h"
#include "subVector.h"
#include "subMatrix.h"
#include "matrixTools.h"
#include "cdsAuto_ptr.h"
#include "cdsIomanip.h"

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
using CDSMath::sq;
#endif /* USE_CDS_NAMESPACE */

static CDSMat33 makeIdentityRotation();
static const CDSMat33 R_I = makeIdentityRotation(); // handy to have around


void
InertiaTensor::calc(const CDSVec3&     center,
                    const AtomList& aList) 
{
    set(0.0);
    CDSMat33 &m = *this;
    for (int cnt=0 ; cnt<aList.size() ; cnt++) {
        const IVMAtom* a = &(*aList[cnt]);
        m(0,0) += a->mass * (sq(a->pos(1)-center(1)) + sq(a->pos(2)-center(2)));
        m(1,1) += a->mass * (sq(a->pos(0)-center(0)) + sq(a->pos(2)-center(2)));
        m(2,2) += a->mass * (sq(a->pos(0)-center(0)) + sq(a->pos(1)-center(1)));
        m(0,1) -= a->mass * (a->pos(0)-center(0)) * (a->pos(1)-center(1));
        m(1,2) -= a->mass * (a->pos(1)-center(1)) * (a->pos(2)-center(2));
        m(0,2) -= a->mass * (a->pos(0)-center(0)) * (a->pos(2)-center(2));
    }
    m(1,0) = m(0,1);
    m(2,1) = m(1,2);
    m(2,0) = m(0,2);
}

////////////////////////////////////////////////
// Implementation of AtomClusterNode methods. //
////////////////////////////////////////////////

// Creates a node with just the 'inboard joint' attachment atom
// and links it to its parent.
AtomClusterNode::AtomClusterNode(const IVM*       ivm,
                                 IVMAtom*         hingeAtom,
                                 const IVMAtom*   parAtom,
                                 AtomClusterNode* parNode) 
  : ivm(ivm),
    atoms(0,1),
    parentAtom(parAtom),
    stateOffset(-1),
    jointType(UnknownJointType),
    jointIsReversed(false),
    parent(parNode),
    children(0,0),
    level(hingeAtom->index==0 ? 0 : parNode->getLevel()+1),
    rbIndex(-1)
{
    atoms.append(hingeAtom);
    hingeAtom->node = this;
}

void AtomClusterNode::addChild(AtomClusterNode* node) {
    children.append( node );
}

void AtomClusterNode::calcAtomPos(const CDSMat33& R_GB, const CDSVec3& OB_G) {
    atoms[0]->pos = OB_G;
    for (int i=1 ; i<atoms.size() ; i++) {
        IVMAtom& a = *atoms[i];
        a.station_G = R_GB * a.station_B;
        a.pos       = OB_G + a.station_G;
    }
}

void AtomClusterNode::calcAtomVel(const Vec6& V_OB_G) {
    const CDSVec3& w_G = *reinterpret_cast<const CDSVec3*>(&V_OB_G[0]);
    const CDSVec3& v_G = *reinterpret_cast<const CDSVec3*>(&V_OB_G[3]);

    atoms[0]->vel = v_G;
    for (int i=1 ; i<atoms.size() ; i++) {
        IVMAtom& a = *atoms[i];
        a.vel = v_G + cross(w_G, a.station_G);
    }
}

void AtomClusterNode::calcSpatialForce(Vec6& F_OB_G) {
    CDSVec3 moment(0.), force(0.);
    // notice that the sign is screwey [??? CDS comment, I don't see it (sherm)]
    for (int i=0 ; i<atoms.size() ; i++) {
        const IVMAtom& a = *atoms[i];
        CDSVec3 aForce = a.force;
        if (ivm->frictionCoeff() != 0.)
            aForce += a.fric * a.vel * (1. - ivm->bathTemp()/ivm->currentTemp());
        moment += cross(a.pos - atoms[0]->pos, aForce);
        force  += aForce;
    }
    F_OB_G = blockVec(moment, force);
}

// Given a set of desired atomic velocities, combine these into a mass-weighted
// spatial impulse for this node. This can then be used in the dynamic equations
// instead of a force, yielding velocities rather than accelerations.
void AtomClusterNode::calcSpatialImpulse(Vec6& Impulse_OB_G) {
    CDSVec3 angular(0.), linear(0.);
    for (int i=0 ; i<atoms.size() ; i++) {
        const IVMAtom& a = *atoms[i];
        const CDSVec3 aMomentum = a.mass*a.vel;
        angular += cross(a.pos - atoms[0]->pos, aMomentum);
        linear  += aMomentum;
    }
    Impulse_OB_G = blockVec(angular, linear);
}

ostream& 
operator<<(ostream& s, const AtomClusterNode& node) {
    for (int i=0 ; i<node.atoms.size() ; i++) {
        IVMAtom* a = node.atoms[i];
        s << a << "p=" << a->pos << "v=" << a->vel << "| ";
    }
    s << "(rigid body " << node.getRBIndex() << ")";
    return s;
}

//////////////////////////////////////////////////
// Define classes derived from AtomClusterNode. //
//////////////////////////////////////////////////

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class GroundBody : public AtomClusterNode {
public: 
    GroundBody(const AtomClusterNode* node)
      : AtomClusterNode(*node)
    { 
        jointType = ThisIsGround;
        jointIsReversed = false;

        for (int i=0 ; i<atoms.size() ; i++) atoms[i]->vel = CDSVec3(0.0); 
        for (l_int i=0 ; i<atoms.size() ; i++) atoms[i]->node = this; 
    }
    ~GroundBody() {}

    // Implementations of virtual AtomClusterNode methods.
    const char* type() const { return "ground"; }
    void print(int) {}
};

/**
 * This is the templatized abstract base class from which all the joint classes are
 * derived.
 */
template<int dof>
class AtomClusterNodeSpec : public AtomClusterNode {
public:

    // We're presented with a node (Body+Joint) in the reference configuration. Derive appropriate
    // quantities for the node, such as the atom stations on the body, mass properties
    // in the body frame, and the reference locations and orientations.
    //
    // XXX Ball and free joints require that we insert a dummy 'atom' as the origin if
    //     this is a base node and it is not bonded to ground (unless there is already
    //     a dummy atom as atom[0]). This is chosen to be
    //     a point which is the geometric center of the entire outboard tree, but
    //     attached to the current body. At the moment I'm just preserving the
    //     existing behavior; later I hope to trash it. (sherm)
    //
    AtomClusterNodeSpec(const AtomClusterNode* node, int& cnt, const CDSMat33& rotBJ,
                        JointType jt, bool isReversed,
                        bool addDummyOrigin=false)
      : AtomClusterNode(*node)
    {
        stateOffset = cnt;
        cnt+=dof;   // leave room for this node's state variables

        jointType       = jt;
        jointIsReversed = isReversed;


        for (int i=0 ; i<atoms.size() ; i++)
            atoms[i]->node = this; 

        // If requested, check whether we need a dummy origin "atom" and if so
        // insert one at the beginning of the atoms list.
        if (addDummyOrigin)
            addDummyOriginIfNeeded();

        // We *may* have a new atoms[0] now.

        refBinP.setFrame(R_I, atoms[0]->pos - getParent()->getAtom(0)->pos);
        JinB.setFrame(rotBJ, CDSVec3(0.)); // joint is always located at body origin
        calcBodyProperties();
    }

    virtual ~AtomClusterNodeSpec() {}


    /*virtual*/int  getDOF() const { return dof; }
    virtual int  getDim() const { return dof; } // dim can be larger than dof
    virtual void print(int);

private:
    void addDummyOriginIfNeeded();
    void calcBodyProperties();

    // This serves as the owner for an extra 'dummy' IVMAtom needed by some joints.
    CDS::auto_ptr<IVMAtom> cmAtom;

};


//////////////////////////////////////////
// Derived classes for each joint type. //
//////////////////////////////////////////

static CDSMat33 makeJointFrameFromZAxis(const CDSVec3& zVec);


/**
 * Translate (Cartesian) joint. This provides three degrees of translational freedom
 * which is suitable (e.g.) for connecting a free atom to ground.
 * The joint frame J is aligned with the body frame B.
 */
class ACNodeTranslate : public AtomClusterNodeSpec<3> {
public:
    /*virtual*/ const char* type() const { return "translate"; }

    ACNodeTranslate(const AtomClusterNode* node, int& cnt)
      : AtomClusterNodeSpec<3>(node,cnt,R_I,CartesianJoint,false)
    { }
};

/**
 * Free joint. This is a six degree of freedom joint providing unrestricted 
 * translation and rotation for a free rigid body.
 * The joint frame J is aligned with the body frame B.
 */
class ACNodeTranslateRotate3 : public AtomClusterNodeSpec<6> {
public:
    ACNodeTranslateRotate3(const AtomClusterNode* node,
                          int&             cnt,
                          bool             shouldUseEuler)
      : AtomClusterNodeSpec<6>(node,cnt,R_I,FreeJoint,false,true), 
        useEuler(shouldUseEuler)
    { 
        if ( !useEuler )
            cnt++;
    }

    /*virtual*/ const char* type() const { return "full"; }
    /*virtual*/ int  getDim() const { return useEuler ? 6 : 7; } 
private:
    bool   useEuler;          //if False, use Quaternion rep.
};

/**
 * Ball joint. This provides three degrees of rotational freedom, i.e.,
 * unrestricted orientation.
 * The joint frame J is aligned with the body frame B.
 */
class ACNodeRotate3 : public AtomClusterNodeSpec<3> {
public:
    ACNodeRotate3(const AtomClusterNode* node,
                  int&                   cnt,
                  bool                   shouldUseEuler)
      : AtomClusterNodeSpec<3>(node,cnt,R_I,OrientationJoint,false,true), 
        useEuler(shouldUseEuler)
    { 
        if ( !useEuler )
            cnt++;
    }

    /*virtual*/ const char* type() const { return "rotate3"; }
    /*virtual*/ int  getDim() const { return useEuler ? 3 : 4; }
private:
    bool   useEuler;          //if False, use Quaternion rep.
};

/**
 * U-joint like joint type which allows rotation about the two axes
 * perpendicular to zDir. This is appropriate for diatoms and for allowing 
 * torsion+bond angle bending.
 */
class ACNodeRotate2 : public AtomClusterNodeSpec<2> {
public:
    ACNodeRotate2(const AtomClusterNode* node,
                 const CDSVec3&             zVec,
                 int&                    cnt)
      : AtomClusterNodeSpec<2>(node,cnt,makeJointFrameFromZAxis(zVec),
                               UJoint,false)
    { 
    }

    /*virtual*/ const char* type() const { return "rotate2"; }
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis.
 */
class ACNodeTranslateRotate2 : public AtomClusterNodeSpec<5> {
public:
    ACNodeTranslateRotate2(const AtomClusterNode*  node,
                          const CDSVec3&              zVec,
                          int&                     cnt)
      : AtomClusterNodeSpec<5>(node,cnt,makeJointFrameFromZAxis(zVec),
                               FreeLineJoint,false)
    { 
    }

    /*virtual*/ const char* type() const { return "diatom"; }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis.
 */
class ACNodeTorsion : public AtomClusterNodeSpec<1> {
public:
    ACNodeTorsion(const AtomClusterNode*   node,
                 const CDSVec3&               rotDir,
                 int&                      cnt)
      : AtomClusterNodeSpec<1>(node,cnt,makeJointFrameFromZAxis(rotDir),
                               TorsionJoint,false)
    { 
    }

    /*virtual*/ const char* type() const { return "torsion"; }
};


//void 
//combineNodes(const AtomClusterNode* node1,
//	       const AtomClusterNode* node2)
//  //
//  // group the atoms of node2 with those of node 1
//  // 
//  //
//{
// cout << "combineNodes: merging node containing atoms: { ";
// cout << node1->atoms[0]->index;
// for (int i=0 ; i<node1->atoms.size() ; i++)
//   cout << " , " << node1->atoms[i]->index ;
// cout << "}\n";
// cout << "\twith that containing atoms:";
// cout << node2->atoms[0]->index;
// for (l_int i=0 ; i<node2->atoms.size() ; i++)
//   cout << " , " << node2->atoms[i]->index ;
// cout << "}\n";
//
// int groupExists=0;
// Atom* remAtom = node1->atoms[0];
// using InternalDynamics::groupList;
// for (l_int i=0 ; i<groupList.size() ; i++)
//   if ( groupList[i].contains(remAtom->index) ) {
//     groupExists=1;
//     for (int j=0 ; j<node2->atoms.size() ; j++)
//	 groupList[i].append( node2->atoms[j]->index );
//     break;
//   }
//
// if ( !groupExists ) {
//   CDSList< int > l;
//   l.append( remAtom->index );
//   for (int j=0 ; j<node2->atoms.size() ; j++)
//     l.append( node2->atoms[j]->index );
//   groupList.append( l );
// }
//} /* combineNodes */

//
// Given a passed-in proto-AtomClusterNode, generate one that has a joint
// as close as possible to the requested one. On entry, cnt indicates
// the next available state variable slot; we grab some and update cnt
// on exit.
//
/*static*/ AtomClusterNode*
AtomClusterNode::constructFromPrototype
    (AtomClusterNode*&                  node,
     const InternalDynamics::HingeSpec& hingeSpec,
     int&                               cnt)
{
    using InternalDynamics::Exception;

    AtomClusterNode* newNode=0;
    const IVM* ivm = node->ivm;
    CDSString type = hingeSpec.type;
    type.downcase();

    if ( type == "ground" ) 
        newNode = new GroundBody(node);

    if (   (type == "bendtorsion" && node->level>2)
        || (type == "bend" && !node->isGroundNode())
              && (node->atoms.size()>1 || node->children.size()>0)
              &&  node->atoms[0]->bonds.size()>0 ) 
    {
        IVMAtom* atom0 = ivm->getAtoms()[ hingeSpec.atom0 ];
        IVMAtom* atom1 = ivm->getAtoms()[ hingeSpec.atom1 ];

        //direction perp. to two bonds
        CDSVec3 dir = cross(atom0->pos - node->atoms[0]->pos,
                         atom1->pos - node->atoms[0]->pos);
        if ( norm(dir) > 1e-12 ) {
            newNode = new ACNodeTorsion(node, dir, cnt);
        } else {
            // the two bonds are colinear - this is an error unless
            //   node has exactly two bonds
            //   and one of atom0 or atom1 has only one bond.
            if ( node->atoms[0]->bonds.size()==2
                 && (atom0->bonds.size() == 1 || atom1->bonds.size() == 1))
            {
                //in this case, the plane of the bend angle needs only to be 
                // perpendicular to the bond
                //rotation matrix which takes the z-axis
                // to (atoms[0]->pos - parentAtom->pos)
                // notes of 12/6/99 - CDS
                const CDSVec3 u = unitVec( atom1->pos - node->atoms[0]->pos );
                double theta = acos( u.z() );
                double psi   = atan2( u.x() , u.y() );
                double a[] = { cos(psi) , cos(theta)*sin(psi) , sin(psi)*sin(theta),
                              -sin(psi) , cos(theta)*cos(psi) , cos(psi)*sin(theta),
                               0        , -sin(theta)         , cos(theta)         };
                CDSVec3 x(1,0,0);
                dir = CDSMat33(a) * x;
                cerr << "norm: " << dot(u,dir) << endl; // should be zero
                newNode = new ACNodeTorsion(node, dir, cnt);
            }
        }
    } 

    if ( type == "torsion" ) {
        if (node->isBaseNode() &&             // allow full dof for base node
            node->parentAtom->index==0)       // unless it is bonded to fixed atom
            //FIX: this seems ill-thought out
            type = "full"; 
        else if ( node->atoms.size()>1 || node->children.size()>0 )
            newNode = new ACNodeTorsion(node, 
                                       node->atoms[0]->pos - node->parentAtom->pos,
                                       cnt);
    } 

    if (type == "rotate" ) {
        if ( node->atoms.size()>2 )
            newNode = new ACNodeRotate3(node,cnt,
                                       ivm->minimization());
        else if ( node->atoms.size()==2 )
            newNode = new ACNodeRotate2(node,
                            node->atoms[1]->pos - node->atoms[0]->pos,cnt);
    }

    if (type == "translate")
        newNode = new ACNodeTranslate(node,cnt);

    if ( type == "full" || type == "unknown") { // all other cases- full freedom of movement
        if ( node->atoms.size()>2 || node->atoms.size()==2 && node->children.size() )
            newNode = new ACNodeTranslateRotate3(node,cnt,
                                                ivm->minimization());
        else if ( node->atoms.size()==2 )
            newNode = new ACNodeTranslateRotate2(node,
                            node->atoms[1]->pos - node->atoms[0]->pos, cnt);
        else if ( node->atoms.size()==1 )
            newNode = new ACNodeTranslate(node,cnt);
    }

    if ( !newNode ) {
        cerr << "Bad Hinge type or topology not supported.\n";
        cerr << "atom cluster node at level " << node->level 
             << " containing atoms: " << node->atoms << '\n';
        throw Exception("Bad Hinge type or topology not supported");
    }

    //
    // connect to other related hinge
    //
    if (!node->isGroundNode()) {
        int index = node->parent->children.getIndex(node);
        node->parent->children[index] =  newNode;
    }
    for (int i=0 ; i<node->children.size() ; i++)
        node->children[i]->parent = newNode;

    delete node; // now we're done with this.
    return newNode;
}


///////////////////////////////////////////////////////////////
// Implementation of AtomClusterNodeSpec base class methods. //
///////////////////////////////////////////////////////////////

// Add a new origin "atom" (station) if certain conditions are met:
//   (1) there is not already a dummy origin, and
//   (2) this is a base (level 1) node, and
//   (3) there is no bond to an atom fixed on ground.
//
// This is a private routine optionally called
// from the AtomClusterNodeSpec constructor.
template<int dof> void 
AtomClusterNodeSpec<dof>::addDummyOriginIfNeeded() {
    if (isBaseNode() &&  atoms[0]->mass != 0.0) {
        bool bound = false;
        for (int i=0 ; i<atoms[0]->bonds.size() ; i++)
            if (atoms[0]->bonds[i]->node->isGroundNode())
                bound = true;

        if ( !bound ) {
            IVMAtom* a = new IVMAtom(-1,0.0);
            cmAtom.reset(a);
            a->pos = AtomTree::findCM( this );
            a->node = this;
            atoms.prepend(a);
        }
    }
}

//
// Calc atom stations, mass, comStation, inertia in B frame (about B origin).
// This is a private routine called from the AtomClusterNodeSpec constructor.
//
template<int dof> void 
AtomClusterNodeSpec<dof>::calcBodyProperties() {
    double  mass = atoms[0]->mass;
    CDSVec3    comStation_B(0.);
    Inertia inertia_OB_B;   // defaults to zero

    const CDSVec3 OB = atoms[0]->pos;
    for (int i=1; i<atoms.size(); ++i) {
        const double   m = atoms[i]->mass;
        const CDSVec3     S = atoms[i]->pos - OB; // atom station

        atoms[i]->station_B  = S;
        mass                += m;
        comStation_B        += m*S;
        inertia_OB_B        += Inertia(m,S);
    }
    comStation_B /= mass;

    massProps.setMassProperties(mass,comStation_B,inertia_OB_B);
}

template<int dof> void
AtomClusterNodeSpec<dof>::print(int verbose) {
    if (verbose&InternalDynamics::printNodeForce) 
        cout << setprecision(8)
             << atoms[0] << ": force: " << atoms[0]->force << '\n';
    if (verbose&InternalDynamics::printNodePos) 
        cout << setprecision(8)
             << atoms[0] << ": pos: " << atoms[0]->pos << ' ' << atoms[0]->vel
             << '\n';
}


/////////////////////////////////////
// Miscellaneous utility routines. //
/////////////////////////////////////

// Calculate a rotation matrix R_BJ which defines the J
// frame by taking the B frame z axis into alignment 
// with the passed-in zDir vector. This is not unique.
// notes of 12/6/99 - CDS
static CDSMat33
makeJointFrameFromZAxis(const CDSVec3& zVec) {
    const CDSVec3 zDir = unitVec(zVec);

    // Calculate spherical coordinates.
    double theta = acos( zDir.z() );             // zenith (90-elevation)
    double psi   = atan2( zDir.x() , zDir.y() ); // 90-azimuth

    // This is a space fixed 1-2-3 sequence with angles
    // a1=-theta, a2=0, a3=-psi. That is, to get from B to J
    // first rotate by -theta around the B frame x axis, 
    // then rotate by -psi around the B frame z axis. (sherm)

    const double R_BJ[] = 
        { cos(psi) , cos(theta)*sin(psi) , sin(psi)*sin(theta),
         -sin(psi) , cos(theta)*cos(psi) , cos(psi)*sin(theta),
          0        , -sin(theta)         , cos(theta)         };
    return CDSMat33(R_BJ); // == R_PJi
}

static CDSMat33
makeIdentityRotation() {
    CDSMat33 ret(0.);
    ret.setDiag(1.);
    return ret;
}
