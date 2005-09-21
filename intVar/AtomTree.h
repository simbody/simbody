#ifndef ATOM_TREE_H_
#define ATOM_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

#include "RigidBodyTree.h"

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
 *
 * TODO Currently constraints are dealt with exclusively at the AtomTree level (using RigidBodyTree
 * services) because the constraint handling has not yet been moved into the RigidBodyTree
 * and is still entangled with atoms. (sherm)
 */
class AtomTree {
public:
    IVM*                         ivm;       // owner of the atoms
    CDSList<AtomClusterNodeList> nodeTree;  // the atom cluster tree
    CDSList<AtomLoop>            loops;

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
    int getNRigidBodies() const { return rbTree.getNBodies(); }

    // Dig up cluster info from rigid body tree.

    double getClusterMass(int level, int indx) const;

    /// Get cluster's center of mass, measured and expressed in the the
    /// ground frame. Requires previous call to setPos().
    const Vec3& getClusterCOM_G(int level, int indx) const;

    /// Get the spatial velocity of this cluster. Requires previous call to setVel().
    const Vec6& getClusterSpatialVel(int level, int indx) const;

    /// Set this cluster's inboard joint coordinates to best approximate
    /// the desired spatial velocity, taking into account the spatial
    /// velocity of the parent.
    void setClusterVelFromSVel(int level, int indx, const Vec6& sVel);

    /// Calculate the kinetic energy contribution of a single cluster
    /// from its spatial velocity. Requires previous call to setVel().
    double calcClusterKineticEnergy(int level, int indx) const;

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos) { rbTree.setPos(pos); calcAtomPos(); }
    void setVel(const RVec& vel) { rbTree.setVel(vel); calcAtomVel(); }
    void setPosVel(const RVec& p, const RVec& v) { setPos(p); setVel(v); }

    RVec getPos() const { RVec pos(getIVMDim()); rbTree.getPos(pos); return pos; }
    RVec getVel() const { RVec vel(getIVMDim()); rbTree.getVel(vel); return vel; }

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
    void calcP() { rbTree.calcP(); }

    /// Dynamics -- calculate articulated body remainder forces.
    void calcZ() {
        calcSpatialForces();
        rbTree.calcZ(spatialForces);
    }

    void applyForces(const VecVec6& spatialForces) {
        rbTree.calcZ(spatialForces);
    }
    
    /// Use dynamics equations to solve velocity problem rather than accelerations.
    /// Normal use requires that setVel(0) has been done earlier, but we won't check.
    void propagateSVel() {
        calcSpatialImpulses(); // sets spatialForces vector
        rbTree.calcZ(spatialForces);
    }

    void calcPandZ() { calcP(); calcZ(); }
    void calcY() { rbTree.calcY(); }

    /// Recalculate unconstrained accelerations given a new set of forces
    /// at the same state.
    void updateAccel() { 
        calcSpatialForces();
        rbTree.calcTreeForwardDynamics(spatialForces);
    }

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
    void calcSpatialImpulses();
    int  getIVMDim() const;  // pull 'dim' from ivm
};

#endif /* ATOM_TREE_H_ */
