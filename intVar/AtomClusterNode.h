#ifndef ATOM_CLUSTER_NODE_H_
#define ATOM_CLUSTER_NODE_H_

#include "internalDynamics.h"

#include "vec3.h"
#include "Mat33.h"
#include "cdsList.h"

class IVM;
class IVMAtom;
class AtomClusterNode;
class RigidBodyNode;
typedef CDSList<IVMAtom*>         AtomList;
typedef CDSList<AtomClusterNode*> AtomClusterNodeList;
typedef FixedVector<double,6>     Vec6;

class InertiaTensor : public Mat33 {
public:
    InertiaTensor() : Mat33(0.0) {}
    //  InertiaTensor(const InertiaTensor&);
    void calc(const Vec3&     center,
              const AtomList&       );
};

/**
 * First crack at separating model building from execution (sherm).
 *
 * Every AtomClusterNode contains a list of pointers to the atoms rigidly affixed to the
 * body associated with that node. The 0th atom in that list is the one which is
 * connected to the node's parent by a joint, to a specific atom on the parent
 * node stored as 'parentAtom' here. Our 0th atom's body-fixed station serves as
 * the origin for the body frame, and thus has coordinates [0,0,0] in that frame.
 */
class AtomClusterNode {
public:
    class VirtualBaseMethod {};    // an exception

    AtomClusterNode(const IVM*        ivm,
                    IVMAtom*          hingeAtom,
                    const IVMAtom*    parentAtom,
                    AtomClusterNode*  parentNode);
    AtomClusterNode(const AtomClusterNode&);
    AtomClusterNode& operator=(const AtomClusterNode&);
    virtual ~AtomClusterNode() {}

    void addChild(AtomClusterNode*);

    const AtomClusterNode* getParent()     const {return parent;}
    AtomClusterNode*       updParent()           {return parent;}
    const IVMAtom*         getParentAtom() const {return parentAtom;}

    int            getNAtoms()     const {return atoms.size();}    
    const IVMAtom* getAtom(int i)  const {return (i<atoms.size()?atoms[i]:0); }

    int                    getNChildren()  const {return children.size();}
    const AtomClusterNode* getChild(int i) const {return (i<children.size()?children[i]:0);}
    AtomClusterNode*       updChild(int i)       {return (i<children.size()?children[i]:0);}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the Ground node at level 0. Level 1 nodes (directly connected to the
    /// Ground node) are called 'base' nodes.
    int              getLevel()  const {return level;}
    bool             isGroundNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }

    int              getStateOffset() const {return stateOffset;}

    int  getRBIndex()       const { return rbIndex; }
    void setRBIndex(int ix)       { rbIndex=ix; }
    
    /// From a temporary node which has been used to collect up clusters of atoms,
    /// generate one that includes a joint, and free the old one.
    static AtomClusterNode* constructFromPrototype
        (AtomClusterNode*&                  oldNode,
         const InternalDynamics::HingeSpec& type,
         int&                               cnt);

    /// Given a spatial orientation and location for this cluster, calculate
    /// where all the atoms are in space.
    void calcAtomPos(const Mat33& R_GB, const Vec3& OB_G);

    /// Given the spatial velocity (angular,linear) of this cluster at its origin,
    /// calculate the velocity of each atom and store that with the atom. 
    /// calcAtomPos() must have been called previously.
    void calcAtomVel(const Vec6& V_OB_G);

    /// Assuming calcAtomPos() and calcAtomVel() have already been called, and
    /// that the atoms have forces applied to them, this will combine all those
    /// atomic forces into a single spatial force acting on the cluster at its
    /// origin and return that force. TODO: currently the bath temperature maintenance force is
    /// generated here and added to the total.
    void calcSpatialForce(Vec6& F_OB_G);

    /// Given a set of desired atomic velocities, combine these into a mass-weighted
    /// spatial impulse for this node. This can then be used in the dynamic equations
    /// instead of a force, yielding velocities rather than accelerations.
    void calcSpatialImpulse(Vec6& Impulse_OB_G);

    virtual int getDOF() const {return 0;} //number of independent dofs
    virtual int getDim() const {return 0;} //# of generalized coords (>=#dofs)

    virtual const char* type() { return "unknown"; }
    virtual void print(int) { throw VirtualBaseMethod(); }
public:
    const IVM*          ivm;
    AtomList            atoms;

private:
    const IVMAtom*      parentAtom;   // atom in parent to which hinge is attached

    int                 stateOffset;
    AtomClusterNode*    parent; 
    AtomClusterNodeList children;
    int                 level;        // how far from base

    // This maps this atom cluster to a rigid body node owned by
    // the same AtomTree that owns the cluster.
    int                 rbIndex;

    friend void combineNodes(const AtomClusterNode* node1,
                             const AtomClusterNode* node2);

    //  static void groupTorsion(const HingeNode*);

    friend ostream& operator<<(ostream& s, const AtomClusterNode&);

    friend class AtomTree;
    friend class AT_Build;
};

#endif /* ATOM_CLUSTER_NODE_H_ */
