#ifndef ATOM_TREE_H_
#define ATOM_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

class Vec3;
class HingeNode;
typedef CDSList<HingeNode*>   HingeNodeList;
typedef CDSVector<double,1>   RVec;   // first element has index 1

/**
 * The RigidBodyTree class owns the tree of joint-connected rigid bodies, called
 * HingeNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * RigidBodyTree is the owner of the HingeNode objects (which are abstract), pointers to
 * which are stored in the nodeTree.
 */
class RigidBodyTree {
public:
    RigidBodyTree() { }

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is not necessarily a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addNode(HingeNode* node) ;
    /*{
        const int level = node->level;
        if (nodeTree.size()<=level) nodeTree.resize(level+1);
        const int offset = nodeTree[level].size();
        nodeTree[level].append(node);
        const int nodeNum = nodeNum2NodeMap.size();
        nodeNum2NodeMap.append(HingeNodeIndex(level,offset));
        node->setNodeNum(nodeNum);
        return nodeNum;
    }*/

    // includes ground
    int getNBodies() const { return nodeNum2NodeMap.size(); }

    int getLevel(int nodeNum) const { return getHingeNode(nodeNum)->getLevel(); }
    
    void velFromCartesian(const RVec& pos,
                                RVec& vel);
    void calcP();
    void calcZ(); 
    void calcPandZ();

    int getDOF(); 
    int getDim(); 
    // deallocate given molecule
    void destructNode(HingeNode*); 

    void setPosVel(const RVec& pos,
                   const RVec& vel);
    void setVel(const RVec& vel);
    void enforceConstraints(RVec& pos,
                            RVec& vel);

    RVec getPos() const;
    RVec getVel() const;
    RVec calcGetAccel();
    RVec getAccel();
    void updateAccel();
    RVec getInternalForce();
    void propagateSVel();
    
    static void addCM(const HingeNode* n,
                      double&          mass,
                      Vec3&            pos);

    static Vec3 findCM(const HingeNode* n);

    const HingeNode& getHingeNode(int nodeNum) const {
        const HingeNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *nodeLevels[ix.level][ix.offset];
    }
    HingeNode& updHingeNode(int nodeNum) {
        const HingeNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *nodeLevels[ix.level][ix.offset];
    }

private:
    struct HingeNodeIndex {
        HingeNodeIndex(int l, int o) : level(l), offset(o) { }
        int level, offset;
    };

    // This holds nodes and serves to map (level,offset) to nodeSeqNo.
    CDSList<HingeNodeList> nodeLevels;

    // Map nodeNum to (level,offset).
    CDSList<HingeNodeIndex> nodeNum2NodeMap;
};


class AtomClusterNode;
class IVM;
class IVMAtom;
class AT_Build;
typedef CDSList<IVMAtom*>     AtomList;
typedef CDSList<AtomClusterNode*>   AtomClusterNodeList;

/**
 * AtomTree owns the RigidBodyTree and maintains a partitioning of IVM's atoms onto the
 * rigid bodies.
 */
class AtomTree {
    IVM*                ivm;
    RigidBodyTree       rbtree;
    CDSList<AtomClusterNodeList> nodeTree;
public:
    AtomTree(IVM*);
    ~AtomTree();
    AtomTree(const AtomTree&);
    AtomTree& operator=(const AtomTree&);

    void addAtomClusterNode(AtomClusterNode* node);

    void addMolecule(AtomClusterNode* groundNode,
                     IVMAtom*         atom      ) ;
    void markAtoms(CDSVector<bool,0>& assignedAtoms);

    friend ostream& operator<<(ostream&,const AtomTree&);
};

#endif /* ATOM_TREE_H_ */
