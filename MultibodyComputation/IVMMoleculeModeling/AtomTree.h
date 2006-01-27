#ifndef ATOM_TREE_H_
#define ATOM_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

#include "IVMRigidBodyTree.h"

#include "IVMMoleculeRBTreeInterface.h"

class AtomClusterNode;
class LengthConstraints;
class IVM;
class IVMAtom;
class AT_Build;

typedef CDSList<IVMAtom*>         AtomList;
typedef CDSList<AtomClusterNode*> AtomClusterNodeList;


class AtomLoop {
public:
    AtomLoop()
      : tip1(0), tip2(0), rbDistConstraintIndex(-1) {}
    AtomLoop(IVMAtom* baseAtom, IVMAtom* tipAtom)
      : tip1(baseAtom), tip2(tipAtom), rbDistConstraintIndex(-1) { }

    IVMAtom* getTip1() const { assert(isValid()); return tip1; }
    IVMAtom* getTip2() const { assert(isValid()); return tip2; }

    void setRBDistanceConstraintIndex(int ix) { rbDistConstraintIndex=ix; }
    int  getRBDistanceConstraintIndex() const {
        assert(isValid() && rbDistConstraintIndex >= 0);
        return rbDistConstraintIndex;
    }
    bool isValid() const { return tip1 && tip2; }

private:
    IVMAtom* tip1;
    IVMAtom* tip2;
    int      rbDistConstraintIndex; // remember the corresponding constraint
};

/**
 * AtomTree owns the atom cluster tree, partioning the IVM-owned atoms into interconnected
 * clusters, and a RigidBodyTree in which each rigid body corresponds to one
 * of the clusters. AtomTree performs useful computations by delegation to
 * the RigidBodyTree and mapping the results back on to the clusters.
 */
class AtomTree {
public:
    IVM*                         ivm;       // owner of the atoms
    CDSList<AtomClusterNodeList> nodeTree;  // the atom cluster tree
    CDSList<AtomLoop>            loops;

private:
    IVMMoleculeRBTreeInterface*  rbTree;     // the pure rigid body tree
    CDSVecVec6                   spatialForces;

    const IVMMoleculeRBTreeInterface& getRBTree() const 
      { assert(rbTree); return *rbTree;}
    IVMMoleculeRBTreeInterface& updRBTree()
      { assert(rbTree); return *rbTree;}

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
    int getNRigidBodies() const { return getRBTree().getNBodies(); }

    // Dig up cluster info from rigid body tree.

    double getClusterMass(int level, int indx) const;

    /// Get cluster's center of mass, measured and expressed in the the
    /// ground frame. Requires previous call to setPos().
    const CDSVec3& getClusterCOM_G(int level, int indx) const;

    /// Get the spatial velocity of this cluster. Requires previous call to setVel().
    const CDSVec6& getClusterSpatialVel(int level, int indx) const;

    /// Set this cluster's inboard joint coordinates to best approximate
    /// the desired spatial velocity, taking into account the spatial
    /// velocity of the parent.
    void setClusterVelFromSVel(int level, int indx, const CDSVec6& sVel);

    /// Calculate the kinetic energy contribution of a single cluster
    /// from its spatial velocity. Requires previous call to setVel().
    double calcClusterKineticEnergy(int level, int indx) const;

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos) { updRBTree().setPos(pos); calcAtomPos(); }
    void setVel(const RVec& vel) { updRBTree().setVel(vel); calcAtomVel(); }
    void setPosVel(const RVec& p, const RVec& v) { setPos(p); setVel(v); }

    RVec getPos() const { RVec pos(getIVMDim()); getRBTree().getPos(pos); return pos; }
    RVec getVel() const { RVec vel(getIVMDim()); getRBTree().getVel(vel); return vel; }

    RVec calcGetAccel();
    RVec getAccel();
    RVec getAccelIgnoringConstraints();
    RVec getInternalForce();

    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const RVec& pos, RVec& vel);

    void fixVel0(RVec& vel); // TODO - get rid of this

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints.
    void enforceConstraints(RVec& pos, RVec& vel);

    /// Dynamics -- calculate articulated body inertias.
    void calcP() { updRBTree().calcP(); }

    /// Dynamics -- calculate articulated body remainder forces.
    void calcZ() {
        calcSpatialForces();
        updRBTree().calcZ(spatialForces);
    }

    void applyForces(const CDSVecVec6& spatialForces) {
        updRBTree().calcZ(spatialForces);
    }

    void calcPandZ() { calcP(); calcZ(); }
    void calcY()     { updRBTree().calcY(); }

    /// Recalculate unconstrained accelerations given a new set of forces
    /// at the same state.
    void updateAccel() { 
        calcSpatialForces();
        updRBTree().calcTreeForwardDynamics(spatialForces);
    }

    int getDOF(); 
    int getDim(); 

    void markAtoms(CDSVector<bool,0>& assignedAtoms);
    static void addCM(const AtomClusterNode* n,
                      double&                mass,
                      CDSVec3&               pos);

    static CDSVec3 findCM(const AtomClusterNode* n);

    friend ostream& operator<<(ostream&,const AtomTree&);

private: 
    void calcAtomPos();
    void calcAtomVel();
    void calcSpatialForces();
    void calcSpatialImpulses();
    int  getIVMDim() const;  // pull 'dim' from ivm
};

#endif /* ATOM_TREE_H_ */
