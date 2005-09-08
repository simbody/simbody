#ifndef ATOM_TREE_H_
#define ATOM_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

#include "RigidBodyTree.h"

class AtomClusterNode;
class IVM;
class IVMAtom;
class AT_Build;
typedef CDSList<IVMAtom*>     AtomList;
typedef CDSList<AtomClusterNode*>   AtomClusterNodeList;

/**
 * AtomTree owns the atom cluster tree, partioning the IVM-owned atoms into interconnected
 * clusters, and a RigidBodyTree in which each rigid body corresponds to one
 * of the clusters. AtomTree performs useful computations by delegation to
 * the RigidBodyTree and mapping the results back on to the clusters.
 *
 * TODO Currently constraints are dealt with exclusively at the AtomTree level (using RigidBodyTree
 * services) because the constraint handling has not yet been moved into the RigidBodyTree
 * and is still entangled with atoms. (sherm)
 */
class AtomTree {
public:
    IVM*                         ivm;       // owner of the atoms
    CDSList<AtomClusterNodeList> nodeTree;  // the atom cluster tree
private:
    RigidBodyTree rbTree;                   // the pure rigid body tree
    VecVec6       spatialForces;
public:
    AtomTree(IVM*);
    ~AtomTree();
    AtomTree(const AtomTree&);
    AtomTree& operator=(const AtomTree&);

    void addAtomClusterNode(AtomClusterNode* node);

    void addMolecule(AtomClusterNode* groundNode,
                     IVMAtom*         atom      ) ;

    /// Call this when all AtomClusterNodes have been built.
    void createRigidBodyTree();

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos) { rbTree.setPos(pos); calcAtomPos(); }
    void setVel(const RVec& vel) { rbTree.setVel(vel); calcAtomVel(); }
    void setPosVel(const RVec& p, const RVec& v) { setPos(p); setVel(v); }

    RVec getPos() const { RVec pos(getIVMDim()); rbTree.getPos(pos); return pos; }
    RVec getVel() const { RVec vel(getIVMDim()); rbTree.getVel(vel); return vel; }
    RVec calcGetAccel() { 
        RVec acc(getIVMDim()); rbTree.calcGetAccel(acc);
        return acc;
    }

    RVec getAccel() {
        calcSpatialForces();
        RVec acc(getIVMDim()); rbTree.getAccel(spatialForces, acc);
        return acc;
    }

    RVec getInternalForce() {
        calcSpatialForces();
        RVec T(getIVMDim()); rbTree.getInternalForce(spatialForces,T);
        return T;
    }

    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const RVec& pos, RVec& vel);

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints.
    void enforceConstraints(RVec& pos, RVec& vel) { rbTree.enforceConstraints(pos,vel); }

    /// Dynamics -- calculate articulated body inertias.
    void calcP() { rbTree.calcP(); }

    /// Dynamics -- calculate articulated body remainder forces.
    void calcZ() {
        calcSpatialForces();
        rbTree.calcZ(spatialForces);
    }

    void calcPandZ() { calcP(); calcZ(); }
    void calcY() { rbTree.calcY(); }

    void updateAccel() { calcSpatialForces(); rbTree.updateAccel(spatialForces); }


    void propagateSVel();

    int getDOF(); 
    int getDim(); 

    // deallocate a node and all its children (recursively)
    void destructNode(AtomClusterNode*);
    void markAtoms(CDSVector<bool,0>& assignedAtoms);
    static void addCM(const AtomClusterNode* n,
                      double&                mass,
                      Vec3&                  pos);

    static Vec3 findCM(const AtomClusterNode* n);

    friend ostream& operator<<(ostream&,const AtomTree&);

private: 
    void calcAtomPos();
    void calcAtomVel();
    void calcSpatialForces();
    int  getIVMDim() const;  // pull 'dim' from ivm
};

#endif /* ATOM_TREE_H_ */
