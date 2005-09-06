#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

#include "RigidBodyNode.h"

class Vec3;
typedef CDSList<RigidBodyNode*>   RBNodePtrList;
typedef CDSVector<double,1>       RVec;   // first element has index 1
typedef CDSList<Vec6>             VecVec6;

/**
 * The RigidBodyTree class owns the tree of joint-connected rigid bodies, called
 * RigidBodyNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * RigidBodyTree is the owner of the RigidBodyNode objects (which are abstract), pointers to
 * which are stored in the tree.
 */
class RigidBodyTree {
public:
    RigidBodyTree() { }
    ~RigidBodyTree();

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addRigidBodyNode(RigidBodyNode* node)
    {
        const int level = node->getLevel();
        if (rbNodeLevels.size()<=level) rbNodeLevels.resize(level+1);
        const int offset = rbNodeLevels[level].size();
        rbNodeLevels[level].append(node);
        const int nodeNum = nodeNum2NodeMap.size();
        nodeNum2NodeMap.append(RigidBodyNodeIndex(level,offset));
        node->setNodeNum(nodeNum);
        return nodeNum;
    }

    // deallocate subtree rooted at the indicated node
    void destructNode(RigidBodyNode*); 

    // includes ground
    int getNBodies() const { return nodeNum2NodeMap.size(); }

//    int getLevel(int nodeNum) const { return getRigidBodyNode(nodeNum)->getLevel(); }

    int getDOF() const; 
    int getDim() const; 

    void realizeParameters();
    void realizeConfiguration(const RVec& pos);
    void realizeMotion(const RVec& vel);
    void realizeAcceleration();

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos);
    void setVel(const RVec& vel);
    void setPosVel(const RVec& pos,
                   const RVec& vel);

    void calcConfigurationKinematics();
    void calcMotionKinematics();

    void getPos(RVec& pos) const;
    void getVel(RVec& vel) const;
    void calcGetAccel(RVec& acc);
    void getAccel(const VecVec6& spatialForces, RVec& acc);
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const RVec& pos, RVec& vel);

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints.
    void enforceConstraints(RVec& pos, RVec& vel);

    /// Dynamics -- calculate articulated body inertias.
    void calcP();

    /// Dynamics -- calculate articulated body remainder forces.
    void calcZ(const VecVec6& spatialForces); 
    void calcPandZ(const VecVec6& spatialForces);

    void calcY();


    void updateAccel(const VecVec6& spatialForces);
    void getInternalForce(const VecVec6& spatialForces, RVec& T);
    void propagateSVel();

    const RigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }
    RigidBodyNode& updRigidBodyNode(int nodeNum) {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

private:
    struct RigidBodyNodeIndex {
        RigidBodyNodeIndex(int l, int o) : level(l), offset(o) { }
        int level, offset;
    };

    // This holds pointers to nodes and serves to map (level,offset) to nodeSeqNo.
    CDSList<RBNodePtrList> rbNodeLevels;

    // Map nodeNum to (level,offset).
    CDSList<RigidBodyNodeIndex> nodeNum2NodeMap;
};

#endif /* RIGID_BODY_TREE_H_ */
