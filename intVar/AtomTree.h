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
 * AtomTree owns the RigidBodyTree and maintains a partitioning of IVM's atoms onto the
 * rigid bodies.
 */
class AtomTree {
public:
    IVM*                ivm;
private:
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

    /// Call this when all AtomClusterNodes have been built.
    void createRigidBodyTree();

    void calcAtomPos(/*state*/);
    void calcAtomVel(/*state*/);
    void calcAtomForces();
    void calcAtomAcc(/*state*/);

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
};

#endif /* ATOM_TREE_H_ */
