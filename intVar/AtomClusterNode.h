#ifndef ATOM_CLUSTER_NODE_H_
#define ATOM_CLUSTER_NODE_H_

#include "internalDynamics.h"

#include "vec3.h"
#include "Mat33.h"
#include "cdsList.h"

class IVM;
class IVMAtom;
class AtomClusterNode;
typedef CDSList<IVMAtom*>         AtomList;
typedef CDSList<AtomClusterNode*> AtomClusterNodeList;

class InertiaTensor : public Mat33 {
public:
    InertiaTensor() : Mat33(0.0) {}
    //  InertiaTensor(const InertiaTensor&);
    void calc(const Vec3&     center,
              const AtomList&       );
};

/**
 * First crack at separating model building from execution.
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

    
    /// From a temporary node which has been used to collect up clusters of atoms,
    /// generate one that includes a joint, and free the old one.
    static AtomClusterNode* constructFromPrototype
        (AtomClusterNode*&                  oldNode,
         const InternalDynamics::HingeSpec& type,
         int&                               cnt);

    virtual int getDOF() const {return 0;} //number of independent dofs
    virtual int getDim() const {return 0;} //# of generalized coords (>=#dofs)
    virtual const char* type() { return "unknown"; }
    virtual void print(int) { throw VirtualBaseMethod(); }
public:
    const IVM*          ivm;
    AtomList            atoms;

private:
    const IVMAtom*      parentAtom;   // atom in parent to which hinge is attached
    AtomClusterNode*    parent; 
    AtomClusterNodeList children;
    int                 level;        // how far from base

    friend void combineNodes(const AtomClusterNode* node1,
                             const AtomClusterNode* node2);

    //  static void groupTorsion(const HingeNode*);

    friend ostream& operator<<(ostream& s, const AtomClusterNode&);

    friend class AtomTree;
    friend class AT_Build;
};

#endif /* ATOM_CLUSTER_NODE_H_ */
