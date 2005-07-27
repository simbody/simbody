#ifndef __dint_node_hh__
#define __dint_node_hh__

#include "phiMatrix.h"
#include "internalDynamics.h"

#include <fixedMatrix.h>
#include <cdsMatrix.h>
#include <cdsVector.h>
#include <cdsList.h>

#include <cassert>

class Vec3;
class IVMAtom;
template<class CHAR> class CDSString;
typedef CDSString<char> String;

typedef CDSList< IVMAtom* >    AtomList;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat6;
typedef CDSMatrix<double>     RMat;
typedef CDSVector<double,1>   RVec;

class InertiaTensor : public Mat3 {
public:
    InertiaTensor() : Mat3(0.0) {}
    //  InertiaTensor(const InertiaTensor&);
    void calc(const Vec3&     center,
              const AtomList&       );
};

class IVM;

/**
 * This is an abstract class representing a body and its (generic) inboard joint, that is,
 * the joint connecting it to its parent. Concrete classes are derived from this one to
 * represent each specific type of joint.
 *
 * HingeNodes are linked into a tree structure, beginning with a special 'origin' node
 * defined to be at level 0 and containing only atoms fixed to ground. The HingeNodes
 * are organized into levels as described in the JMR paper. Level 1 nodes (referred to
 * as 'base nodes') are those attached directly to the origin node, level 2's attach 
 * to level 1's, etc. Every node but the origin has exactly one parent node, whose
 * level is always one less than the current node. Any node may have an arbitrary number 
 * of children, for which it is the unique parent, and all of its children have 
 * level one greater than the current node.
 *
 * Every HingeNode contain an IVMAtom list of the atoms rigidly affixed to the
 * body associated with that node. The 0th atom in that list is the one which is
 * connected to the node's parent by a joint, to a specific atom on the parent
 * node stored as 'parentAtom' here. Our 0th atom's body-fixed station serves as
 * the origin for the body frame, and thus has coordinates [0,0,0] in that frame.
 *
 */
class HingeNode {
public:
    const IVM*    ivm;
    const IVMAtom*   parentAtom;   // atom in parent to which hinge is attached
    class VirtualBaseMethod {};    // an exception

    virtual ~HingeNode() {}
    HingeNode(const IVM*        ivm,
              IVMAtom*          hingeAtom=0,
              const IVMAtom*    parentAtom=0,
              HingeNode*        parentNode=0);
    HingeNode& operator=(const HingeNode&);

    HingeNode*       getParent() const {return parent;}

    /// Return R_GB, the rotation (direction cosine) matrix giving the 
    /// orientation of this body's frame B in the ground frame G.
    const Mat3&      getR_GB()   const {return R_GB;}

    /// Return R_GP, the rotation (direction cosine) matrix giving the
    /// orientation of this body's *parent's* body frame (which we'll call
    /// P here) in the ground frame G.
    const Mat3&      getR_GP()   const {assert(parent); return parent->getR_GB();}

    /// Return this node's level, that is, how many ancestors separate it from
    /// the origin node at level 0. Level 1 nodes (directly connected to the
    /// origin node) are called 'base' nodes.
    int              getLevel()  const {return level;}

    bool             isOriginNode() const { return level==0; }
    bool             isBaseNode()   const { return level==1; }


    const Vec6&      getSpatialVel()   const {return sVel;}
    const Vec6&      getSpatialAcc()   const {return sAcc;}

    const PhiMatrix& getPhi()    const {return phi;}
    const Mat6&      getPsiT()   const {return psiT;}
    const Mat6&      getY()      const {return Y;}

    const IVMAtom*   getAtom(int i)   const {return (i<atoms.size()?atoms[i]:0); }
    const HingeNode* getChild(int i) const {return (i<children.size()?children[i]:0);}

    virtual int           offset() const {return 0;}
    virtual const Vec3&   posCM()  const {throw VirtualBaseMethod();}
    virtual const double& mass()   const {throw VirtualBaseMethod();}

    void addChild(HingeNode*);

    virtual void calcP() {throw VirtualBaseMethod();}
    virtual void calcZ() {throw VirtualBaseMethod();}
    virtual void calcPandZ() {throw VirtualBaseMethod();}
    virtual void calcY() {throw VirtualBaseMethod();}
    virtual void calcInternalForce()        {throw VirtualBaseMethod();}
    virtual void prepareVelInternal()       {throw VirtualBaseMethod();}
    virtual void propagateSVel(const Vec6&) {throw VirtualBaseMethod();}

    virtual int getDOF() const {return 0;} //number of independent dofs
    virtual int getDim() const {return 0;} //dofs plus constraints
    virtual double kineticE() { return 0; }
    virtual double approxKE() { return 0; }
    virtual void setPosVel(const RVec&, const RVec&) {throw VirtualBaseMethod();}
    virtual void setVel(const RVec&)    {throw VirtualBaseMethod();}
    virtual void setVelFromSVel(const Vec6&) {throw VirtualBaseMethod();}
    virtual void enforceConstraints(RVec& pos, RVec& vel) {throw VirtualBaseMethod();}
    virtual const char* type() { return "unknown"; }
    virtual void print(int) { throw VirtualBaseMethod(); }
    virtual void getPos(RVec&) {throw VirtualBaseMethod();}
    virtual void getVel(RVec&) {throw VirtualBaseMethod();}
    virtual void getAccel(RVec&) {throw VirtualBaseMethod();}
    virtual void calcAccel() {throw VirtualBaseMethod();}
    virtual void getInternalForce(RVec&) {throw VirtualBaseMethod();}
    virtual RMat getH() {throw VirtualBaseMethod();}

protected:
    typedef CDSList<HingeNode*>   HingeNodeList;

    HingeNode*    parent; 
    HingeNodeList children;
    int           level;        //how far from base
    AtomList      atoms; 

    Vec6      sVel;        //spatial velocity
    Vec6      sAcc;        //spatial acceleration
    PhiMatrix phi;

    Mat6 psiT;
    Mat6 P;
    Vec6 z;
    Mat6 tau;
    Vec6 Gepsilon;

    // Rotation (direction cosine matrix) expressing the body frame B
    // in the ground frame G. That is, if you have a vector vB expressed
    // body frame and want it in ground, use vG = R_GB*vB. 
    Mat3 R_GB;

    Mat6 Y;  //diagonal components of Omega- for loop constraints

    static const double DEG2RAD; //using angles in degrees balances gradient

    virtual void velFromCartesian() {}
    void removeHinge();

    friend ostream& operator<<(ostream& s, const HingeNode&);
    friend HingeNode* construct(HingeNode*                         oldNode,
                                const InternalDynamics::HingeSpec& type,
                                int&                               cnt);
    friend void combineNodes(const HingeNode* node1,
                             const HingeNode* node2);

    //  static void groupTorsion(const HingeNode*);
    template<int dof> friend class HingeNodeSpec;
    friend class AtomTree;
    friend class AT_Build;
};

#endif /* __dint_node_hh__ */
