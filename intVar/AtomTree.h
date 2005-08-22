#ifndef ATOM_TREE_H_
#define ATOM_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"

class IVM;
class IVMAtom;
class AT_Build;
class Vec3;
class HingeNode;

typedef CDSList< IVMAtom* >     AtomList;
typedef CDSVector<double,1>     RVec;   // first element has index 1

/**
 * The AtomTree class owns the tree of joint-connected rigid bodies, called
 * HingeNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * AtomTree is the owner of the HingeNode objects (which are abstract), pointers to
 * which are stored in the nodeTree.
 *
 * XXX currently the AtomTree and HingeNodes are entangled with atoms in a messy
 * XXX way
 *
 */
class AtomTree { 
public:
    IVM* ivm;
    CDSList< CDSList<HingeNode*> > nodeTree; 

    AtomTree(IVM*);
    ~AtomTree();
    AtomTree(IVMAtom*);
    AtomTree(const AtomTree&);
    AtomTree& operator=(const AtomTree&);

    void addNode(HingeNode* node);
    void addMolecule(HingeNode* groundNode,
                     IVMAtom*   atom      ) ;
    void velFromCartesian(const RVec& pos,
                                RVec& vel);
    void calcP();
    void calcZ(); 
    void calcPandZ();
    void markAtoms(CDSVector<bool,0>& assignedAtoms);
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

    friend ostream& operator<<(ostream&,const AtomTree&);
};

#endif /* ATOM_TREE_H_ */
